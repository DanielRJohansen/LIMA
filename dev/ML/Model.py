from typing import Any, Dict
import torch
import torch.nn as nn

# -----------------------------------------------------------------------------
# 3. Model construction
# -----------------------------------------------------------------------------

def build_encoder(params: Dict[str, Any]) -> nn.Module:
    layers: list[nn.Module] = []
    in_features = 9
    for _ in range(params["num_enc_blocks"]):
        layers.append(nn.Linear(in_features, params["enc_hidden_units"]))
        layers.append(getattr(nn, params["enc_activation"])())
        in_features = params["enc_hidden_units"]
    # final projection to 128-dim before quant
    layers.append(nn.Linear(in_features, 128))
    return nn.Sequential(*layers)


def build_bottleneck(params: Dict[str, Any]) -> nn.Module:
    t = params["bottleneck_type"]
    if t == "SignSTE":
        return SignSTE128()
    elif t == "GumbelBernoulli":
        return GumbelBernoulli128()
    elif t == "VQVAE":
        return VQVAE128()
    else:
        raise ValueError(f"Unknown bottleneck type: {t}")

class SignSTE128(nn.Module):
    """
    Binarizes input to {-1, +1} using sign, with straight-through estimator for gradients.
    """
    def __init__(self) -> None:
        super().__init__()
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # forward: sign; backward: identity gradient approx.
        return torch.sign(x)

class GumbelBernoulli128(nn.Module):
    def __init__(self, tau: float = 1.0) -> None:
        super().__init__()
        self.tau = tau

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # relaxed Bernoulli -> {0,1}
        probs = x.sigmoid()
        gumbel = -torch.empty_like(probs).exponential_().log()
        logits = (torch.log(probs) + gumbel) / self.tau
        return torch.sigmoid(logits)  # soft bits

class VQVAE128(nn.Module):
    def __init__(self, num_codes: int = 256) -> None:
        super().__init__()
        self.codebook = nn.Embedding(num_codes, 128)
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # nearest neighbor lookup
        flat = x.view(-1, 128)
        dist = torch.cdist(flat, self.codebook.weight)
        idx = torch.argmin(dist, dim=1)
        quant = self.codebook(idx).view_as(x)
        return quant


def build_decoder(params: Dict[str, Any]) -> nn.Module:
    layers: list[nn.Module] = []
    in_features = 128
    for _ in range(params["num_dec_blocks"]):
        layers.append(nn.Linear(in_features, params["dec_hidden_units"]))
        layers.append(getattr(nn, params["dec_activation"])())
        in_features = params["dec_hidden_units"]
    # project back to 9 floats
    layers.append(nn.Linear(in_features, 9))
    return nn.Sequential(*layers)