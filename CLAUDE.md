# LIMA

## Coding rules

- Never use `std::format`. Use `Lima::Format` from `Format.h` instead; `std::format` costs ~3 s of compile time per translation unit (see the header comment). In `.cu` files avoid both, since nvcc can hit an internal compiler error on `<format>`; build strings with `std::to_string` there.
