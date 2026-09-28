import cv2
import numpy as np
from PIL import Image

# Load the image
image_path = "table1_new.png"
image = cv2.imread(image_path)

# Convert to grayscale
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# Detect edges in the image
edges = cv2.Canny(gray, 50, 150, apertureSize=3)

# Use HoughLinesP to detect lines
lines = cv2.HoughLinesP(edges, 1, np.pi / 180, threshold=100, minLineLength=100, maxLineGap=10)

# Function to draw horizontal lines back after removing vertical line
def draw_horizontal_lines(image, y, color=(0, 0, 0)):
    cv2.line(image, (0, y), (image.shape[1], y), color, 1)

# Draw the detected lines on the image
if lines is not None:
    vertical_lines = []
    for line in lines:
        x1, y1, x2, y2 = line[0]
        if abs(x1 - x2) < 10:  # Vertical line within a small threshold
            vertical_lines.append((x1, y1, x2, y2))

    # Sort vertical lines by their x coordinate
    vertical_lines.sort()

    # Remove the second vertical line
    if len(vertical_lines) > 2:
        x1, y1, x2, y2 = vertical_lines[2]
        cv2.line(image, (x1, y1), (x2, y2), (255, 255, 255), 3)  # Draw white line to "erase" the vertical line

        # Draw horizontal lines back
        horizontal_lines_y = [y1, y2]
        for y in horizontal_lines_y:
            draw_horizontal_lines(image, y, (0, 0, 0))

# Save the modified image
modified_image_path = "table1_modified.png"
cv2.imwrite(modified_image_path, image)

# Display the modified image
modified_image = Image.open(modified_image_path)
modified_image.show()
