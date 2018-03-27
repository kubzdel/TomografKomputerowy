from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import os
from numpy.fft import fftfreq, fft, ifft

from skimage.filters import gaussian

from line_drawer import \
    bresenham, draw_line, calculatePosition, check_borders, calculatePositionSafe, reconstruct_line
import time
from scipy.ndimage.filters import convolve

def plot_image(image):
    plt.subplot(1, 1, 1), plt.imshow(image, cmap='gray')
    plt.show()

def plot_graph(xAxis, yAxis, title=''):
    import matplotlib.pyplot as plt
    plt.plot(xAxis, yAxis)
    plt.ylabel(title)
    plt.show()

def generate_kernel():
    top = -4/(np.pi**2)
    half = kernel_length // 2
    for i in range(0, kernel_length):
        kernel[i] = ( 0 if (i-half)%2==0 else (top/((i-half)**2)) )
    kernel[kernel_length // 2] = 1

def normalize_image_to_one(image):
    min = np.min(image)
    if min < 0:
        #make all values positive
        image += np.abs(min)
    distance = np.max(image)
    return (image / distance)

def normalize_image(image):
    maxV = np.max(image)
    if maxV == 0:
        return image
    return (image / maxV ) * 255

def parallel_ray(main_point_start, main_angle, minor_angle, diameter, center):
    start = calculatePositionSafe(main_angle + minor_angle, diameter, center)
    end = calculatePositionSafe(main_angle - 180 + (-1)*minor_angle, diameter, center)
    return start, end

def inclined_ray(main_angle, minor_angle, diameter, center):
    end_point = calculatePositionSafe(main_angle + minor_angle, diameter, center)
    return (diameter - end_point[0], diameter - end_point[1])

def get_new_image_shape(old_image):
    if old_image.shape[0] != old_image.shape[1]:
        raise Exception('Image is not a square!')
    #size is sqrt(2)*a -> diagonal of a square
    size = (int)(np.sqrt(2) * max(old_image.shape[0], old_image.shape[1]) + 10)
    return (size, size)


#creates quare matrix that is sqrt(2)*size(image)+10 long
def prepare_image(image):
    size = get_new_image_shape(image)
    new_image = np.zeros(size)
    dy = image.shape[0] // 2
    dx = image.shape[1] // 2
    new_center = size[0] // 2
    # put image into new_image
    new_image[new_center - dy: new_center + dy, new_center - dx: new_center + dx] = image
    return new_image


def clear_before_reconstruction(image):
    image_size = len(image)
    freqencies = fftfreq(image_size).reshape(-1, 1)
    omega = 2 * np.pi * freqencies
    fourier_filter = 2 * np.abs(freqencies)
    if use_cosine_in_fourier:
        fourier_filter *= np.cos(omega)
    projection = fft(image, axis=0) * fourier_filter
    return np.real(ifft(projection, axis=0))


def display_status(num, all):
    #os.system('cls')
    p = num * 100 // all
    print("{status}>{spaces}{percent}%".format(status='-' * (p), spaces=(100 - p - 1) * ' ', percent=p))


def radon(oryg_image, image, angle, n_detectors, sinogram_arr, emission_angle, diameter):
    angle_between_rays = emission_angle / (n_detectors - 1)
    angles = np.arange(-emission_angle/2, emission_angle/2 + emission_angle/n_detectors , angle_between_rays)

    radius = diameter // 2 - 1
    center = (radius, radius)
    start = calculatePositionSafe(angle, diameter, center)
    x = 0
    for i in angles:
        if parallel_rays_mode:
            start, end = parallel_ray(start, angle, i*2, diameter, center)
        else:
            end = inclined_ray(angle, i*2, diameter, center)
        line_sum = draw_line(oryg_image, image, start_point=start, end_point=end)
        sinogram_arr[angle, x] = line_sum
        x += 1

def inverse_radon(image, sinogram, diameter, angle, emission_angle, n_detectors):
    angle_between_rays = emission_angle / (n_detectors - 1)
    angles = np.arange(-emission_angle / 2, emission_angle / 2 + emission_angle / n_detectors, angle_between_rays)

    radius = diameter // 2 - 1
    center = (radius, radius)
    start = calculatePositionSafe(angle, diameter, center)

    x = 0
    for i in angles:
        if parallel_rays_mode:
            start, end = parallel_ray(start, angle, i*2, diameter, center)
        else:
            end = inclined_ray(angle, i*2, diameter, center)
        reconstruct_line(sinogram_value=sinogram[angle, x], reconstruction_image=image, start_point=start, end_point=end)
        x+=1

def rootMeanSquareError(input_image,output_image):
    error = 0
    n = 0
    for x in range(len(input_image)):
        for y in range(len(input_image[x])):
            pixelError = pow(input_image[x][y] - output_image[x][y],2)
            n+=1
            error += pixelError
    error = 1/n * error
    error = np.sqrt(error)
    return error

def darken(input_image):
    for x in range(len(input_image)):
        for y in range(len(input_image[x])):
            if input_image[x][y] <0.4:
                input_image[x][y] = 0
            else:
                input_image[x][y]-=0.4


def process(image):
    #new image
    new_image_size = get_new_image_shape(image)

    #keep original image matrix to get original values
    oryg_image = prepare_image(image)
    reconstructed = np.zeros(new_image_size)
    sinogram_arr = np.zeros((radon_angle, n_detectors))

    #create sinogram
    print('Creating sinogram')

    for i in range(0, radon_angle):
        rays_image = np.zeros(new_image_size)
        display_status(i, radon_angle)
        radon(oryg_image, rays_image, i, n_detectors, sinogram_arr, emission_angle, diameter=new_image_size[0])

    print('Reconstructing image')

    #fourier reconstruction - backprojection
    if use_fourier_reconstruction:
        sinogram_arr = clear_before_reconstruction(sinogram_arr)

    #convolution
    if use_convolution_filter:
        for i in range(0, radon_angle):
            sinogram_arr[i, : ] = convolve(sinogram_arr[i, :], kernel)

    #reconstruct image
    for i in range(0, radon_angle):
        display_status(i, radon_angle)
        inverse_radon(reconstructed, sinogram_arr, diameter=new_image_size[0], angle=i, emission_angle=emission_angle,n_detectors=n_detectors)
        reconstructed = normalize_image_to_one(reconstructed)
        oryg_image = normalize_image_to_one(oryg_image)

        print(rootMeanSquareError(oryg_image,reconstructed))
    if use_convolution_in_output:
        reconstructed = convolve(reconstructed, kernel_reconstructed)
    if use_gauss_in_reconstruction:
        reconstructed = gaussian(reconstructed)

    if normalize_output:
        reconstructed = normalize_image_to_one(reconstructed)
        oryg_image = normalize_image_to_one(oryg_image)


    darken(reconstructed)

    return sinogram_arr, reconstructed,rootMeanSquareError(oryg_image, reconstructed)

#parameters
#   emiters detectors number
n_detectors = 100
#   angle between first and last ray
emission_angle = 40
#   rotation
radon_angle = 180

parallel_rays_mode = True
normalize_output = True

#filters parameters
#sinogram convolution
kernel_length = 5
kernel = np.zeros(( kernel_length))
use_convolution_filter = True

#fourier reconstruction
use_fourier_reconstruction = False
use_cosine_in_fourier = True

#output convolution
use_gauss_in_reconstruction = True
use_convolution_in_output = True
kernel_reconstructed = np.ones((9,9))

#file to transform
file = "photo.png"
watroba = "watroba.png"
directory = os.getcwd() + "\\res\\"
ndetectors = [100]

if __name__ == "__main__":
    generate_kernel()
    errors =[]

    for k in ndetectors:
        n_detectors = k
        image = misc.imread('{dir}{file}'.format(dir=directory, file=file), flatten=True).astype('float64')
        sinogram, reconstructed,error = process(image)
        errors.append(error)
        print(error)
   # plt.plot(ndetectors,errors)
    plt.ylabel("Błąd")
    plt.xlabel("Filtr")
    x = np.arange(4)
    filtry=["brak","splot sinogram","splot wyjście"]
    wyniki=[0.238,0.2375,0.109]
    plt.bar(filtry,wyniki)
    plt.show()
    plot_image(sinogram)
    plot_image(reconstructed)


