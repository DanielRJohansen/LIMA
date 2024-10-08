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
#display resize 220 180
puts "Set display resolution to 1920x1080"


set framesPerStep 2

# Set up view and rotation
display resetview
rotate x by -90
set zoom [expr 1 + 0.8 / double($num_frames * $framesPerStep - 1)]
puts "$zoom"
scale by 2
mol modstyle 0 0 CPK

puts "Start rendering frames $num_frames"
for {set i 0} {$i < $num_frames} {incr i} {
    # Go to the current trajectory frame
    animate goto $i

    for {set j 0} {$j < $framesPerStep} {incr j} {
        # Calculate the overall index of the rendered frame
        set frame_index [expr $i * $framesPerStep + $j]

        scale by $zoom

        # Render the current subframe
        render TachyonInternal "./frames/frame_$frame_index.tga"
    }
}
puts "Script finished"
quit