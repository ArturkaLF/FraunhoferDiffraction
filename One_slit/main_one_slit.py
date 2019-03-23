import math
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pylab
from PIL import Image
from matplotlib.widgets import Button, Slider

"""
@author Arthur Shajkhatarov
@version 1.2

Physical process modeling. Fraunhofer diffraction

1) Tuning():
    a) Wave length. Adjust the color of the model
    b) Slit width.
    c) Screen distance.

2) Plotting():
    a)
    b)
    c)
    
3) Modeling():
    a)
    b)
    c)

# Заполнить эти три поля
# Описать каждую часть программы
# Написать документацию на русском

"""


def wavelength_to_rgb(wavelengthColor, gamma=0.8):
    """
    This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).
    """

    wavelengthColor = float(wavelengthColor)
    if 380 <= wavelengthColor <= 440:
        attenuation = 0.3 + 0.7 * (wavelengthColor - 380) / (440 - 380)
        R = ((-(wavelengthColor - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif 440 <= wavelengthColor <= 490:
        R = 0.0
        G = ((wavelengthColor - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif 490 <= wavelengthColor <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelengthColor - 510) / (510 - 490)) ** gamma
    elif 510 <= wavelengthColor <= 580:
        R = ((wavelengthColor - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif 580 <= wavelengthColor <= 645:
        R = 1.0
        G = (-(wavelengthColor - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif 645 <= wavelengthColor <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelengthColor) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    R *= 255
    G *= 255
    B *= 255
    return [int(R), int(G), int(B)]


def single_slit_diffraction_intensity(slit_width, wavelength, screen_distance, X):
    """
    Takes in slit_width, wavelength, screen distance and a numpy array X(an array of distances from the center).
    Outputs an array of normalized intensities corresponding to X.
  	"""
    return ((np.sin((np.pi * slit_width * X) / (wavelength * screen_distance))) / (
            (np.pi * slit_width * X) / (wavelength * screen_distance))) ** 2


# X values
X = np.arange(-0.007, 0.007, 0.00001)

# The size of the image
imgsize = (1000, 1000)
# Create the image
image = Image.new('RGB', imgsize)

# Standard indicators
wavelength = 500 * (10 ** -9)
slit_width = 100 * (10 ** -6)
screen_distance = 50 * (10 ** -2)

# Y values
Y = single_slit_diffraction_intensity(slit_width, wavelength, screen_distance, X)

plot, = plt.plot(X, Y)
plt.xlabel("Distance from center")
plt.ylabel("Intensity")

# Set coordinates of widgets
axis = (plt.axes([0.75, 0.75, 0.14, 0.05]))
axis2 = (plt.axes([0.75, 0.65, 0.14, 0.05]))
axis3 = (plt.axes([0.75, 0.55, 0.14, 0.05]))

# Set widgets (sliders)
wavelength_slider = Slider(axis, 'Wavelength(nm)', 390, 740, valinit=wavelength * 10 ** 9)
slit_width_slider = Slider(axis2, "Slit Width(micrometers)", 10, 1000, valinit=slit_width * 10 ** 6)
screen_distance_slider = Slider(axis3, "Screen Distance(cm)", 10, 100, valinit=screen_distance * 10 ** 2)


def update(val):
    wavelength = wavelength_slider.val * (10 ** -9)
    slit_width = slit_width_slider.val * (10 ** -6)
    screen_distance = screen_distance_slider.val * (10 ** -2)
    Y = single_slit_diffraction_intensity(slit_width, wavelength, screen_distance, X)
    plot.set_ydata(Y)


wavelength_slider.on_changed(update)
slit_width_slider.on_changed(update)
screen_distance_slider.on_changed(update)


def onButtonAddClicked(event):
    global number

    # Update values ​​after clicking of the button
    wavelength_slider.on_changed(update)
    slit_width_slider.on_changed(update)
    screen_distance_slider.on_changed(update)

    wavelength = wavelength_slider.val * (10 ** -9)
    slit_width = slit_width_slider.val * (10 ** -6)
    screen_distance = screen_distance_slider.val * (10 ** -2)

    # Color at the center
    outerColor = [255, 255, 255]
    # Color at the corners
    # Depends on wavelength
    innerColor = wavelength_to_rgb(wavelength_slider.val)

    for y in range(imgsize[1]):
        for x in range(imgsize[0]):

            # Find the distance to the center
            distanceToCenter = math.sqrt((x - imgsize[0] / 2) ** 2 + (y - imgsize[1] / 2) ** 2)
            X = distanceToCenter / 90000

            distanceToCenter = float(((np.sin((np.pi * slit_width * X) / (wavelength * screen_distance))) / (
                    (np.pi * slit_width * X) / (wavelength * screen_distance))) ** 2)

            if distanceToCenter != distanceToCenter:
                distanceToCenter = 1

            # Calculate r, g, and b values
            r = innerColor[0] * distanceToCenter + outerColor[0] * (1 - distanceToCenter)
            g = innerColor[1] * distanceToCenter + outerColor[1] * (1 - distanceToCenter)
            b = innerColor[2] * distanceToCenter + outerColor[2] * (1 - distanceToCenter)

            # Place the pixel
            image.putpixel((x, y), (int(r), int(g), int(b)))

    # Saving option with his number
    """Need to change the path to the path in your directory"""
    image.save(str(number) + ".jpg")
    img = Image.open('/Users/artur/FraunhoferDiffraction/One_slit/' + str(number) + '.jpg')
    img.show()
    number += 1


# Option counter
number = 1

axes_button_add = pylab.axes([0.7, 0.05, 0.25, 0.075])

# Button creation
button_add = Button(axes_button_add, 'Enter')
button_add.on_clicked(onButtonAddClicked)

# Delete past files

"""Need to change the path to the path in your directory"""
files = glob.glob('/Users/artur/FraunhoferDiffraction/One_slit/*')
files.remove('/Users/artur/FraunhoferDiffraction/One_slit/main_one_slit.py')

if files:
    for f in files:
        os.remove(f)

plt.show()
