# Load your structure file
mol new conf.gro
puts "Loaded conf.gro"

# Load trajectory with waitfor all to ensure full loading
set status [mol addfile trajectory.trr type trr waitfor all]
if {$status != 0} {
    puts "Error loading trajectory file"
    exit
}
puts "Added trajectory.trr with waitfor all"

# Check the number of frames after loading
set num_frames [molinfo top get numframes]
puts "Number of frames after reload: $num_frames"

# Ensure frames directory exists
file mkdir ./frames

# Set image resolution to 1920x1080
display resize 1920 1080
puts "Set display resolution to 1920x1080"

# Set up view and rotation
display resetview
rotate x by -90
set zoom [expr 1 + 0.8 / double($num_frames - 1)]
puts "Zoom factor: $zoom"
scale by 2
mol modstyle 0 0 CPK

# Get the list of existing frames in ./frames
set frame_files [glob ./frames/frame_*.tga]
set last_frame_rendered -1
foreach filename $frame_files {
    # Extract the frame number from the filename
    if {[regexp {frame_(\d+)\.tga} $filename match frame_num_str]} {
        set frame_num [int $frame_num_str]
        if {$frame_num > $last_frame_rendered} {
            set last_frame_rendered $frame_num
        }
    }
}
puts "Last frame rendered: $last_frame_rendered"

# Set starting frame to next frame not yet rendered
set start_frame [expr {$last_frame_rendered + 1}]
puts "Starting from frame: $start_frame"

puts "Start rendering frames from $start_frame to $num_frames"
for {set i $start_frame} {$i < $num_frames} {incr i} {
    # Go to the current trajectory frame
    animate goto $i

    # Apply zoom scaling
    scale by $zoom

    # Render the current frame
    render TachyonInternal "./frames/frame_$i.tga"
}

puts "Script finished"
quit
