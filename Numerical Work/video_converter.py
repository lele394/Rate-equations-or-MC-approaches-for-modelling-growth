import h5py
import numpy as np
import cv2

def create_mp4_from_hdf5(hdf5_filename='simulation_frames_uint8.h5', output_filename='simulation_output.mp4', colors=None, padding_color=(0, 0, 0)):
    """
    Create an MP4 video from frames stored in an HDF5 file using OpenCV.

    Parameters:
    - hdf5_filename: str, name of the HDF5 file containing simulation frames.
    - output_filename: str, name of the output MP4 video file.
    - colors: dict, mapping of frame values to RGB tuples.
    - padding_color: tuple, RGB color for padding (default is black).
    """
    # Default colors if not provided
    if colors is None:
        colors = {
            0: (0, 29, 61),      # Black
            1: (131,208,203),    # light teal
            2: (5,250,124),    # Green
            3: (255, 207, 103),    # orange
            11: (57,124,147), # dark teal
            12: (24, 130, 115), # dark green
            13: (226, 102, 54)  # orange
        }


    # For some reason it looks better
    colors = {
        0: (0, 0, 0),      # Black
        1: (255, 0, 0),    # Red
        2: (0, 255, 0),    # Green
        3: (0, 0, 255),    # Blue
        11: (255, 255, 0), # Yellow
        12: (255, 0, 255), # Magenta
        13: (0, 255, 255)  # Cyan
    }

    # Open the HDF5 file and read the frames
    with h5py.File(hdf5_filename, 'r') as f:
        dataset = f['frames']
        num_frames = dataset.shape[0]
        
        # Create a VideoWriter object using OpenCV
        first_frame = dataset[0, :, :]  # Get the first frame to determine size
        height, width = first_frame.shape
        
        # Ensure the video dimensions are divisible by 16
        new_width = (width + 15) // 16 * 16  # Round up to nearest multiple of 16
        new_height = (height + 15) // 16 * 16  # Round up to nearest multiple of 16
        
        fourcc = cv2.VideoWriter_fourcc(*'MJPG')  # Codec for MP4
        out = cv2.VideoWriter(output_filename, fourcc, 60, (new_width, new_height))  # 60 fps

        # Convert each frame to an image and write to video
        for i in range(num_frames):

            #if i%10 != 0: # skip 100 frame
            #    continue

            # if i<1591:
            #     continue

            if i%10 == 0: print(f'{i}/{num_frames}')
            frame = dataset[i, :, :]  # Get the frame
            
            # Create a new RGB image
            img_array = np.zeros((height, width, 3), dtype=np.uint8)  # Initialize an RGB image
            
            # Map values to colors
            for value, color in colors.items():
                # Create a binary mask for the current value
                mask = (frame == value)  # Boolean array where the condition is true
                mask = mask[:, :, np.newaxis]  # Add a new axis for coloring
                
                # Apply the color where the mask is True
                img_array += (mask * np.array(color, dtype=np.uint8))  # Ensure color is uint8
            
            # Ensure the image is padded to the correct size
            padded_img = np.full((new_height, new_width, 3), padding_color, dtype=np.uint8)  # Create a padded image
            padded_img[:height, :width] = img_array  # Place the original image in the top-left corner
            
            out.write(padded_img)  # Write the padded image to the video

    out.release()  # Release the video writer
    print(f"video saved as {output_filename}.")

# Example usage
# create_mp4_from_hdf5('simulation_frames_uint8.h5', 'simulation_output.mp4', padding_color=(255, 255, 255))  # White padding
create_mp4_from_hdf5('cluster_data.h5', 'simulation_output.mp4', padding_color=(0, 0, 0))  # White padding
