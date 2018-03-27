import numpy as np

#returns sum of pixels values that Bresenham's alg went through or reconstructs image along path
def bresenham(start, end, oryg_image = None,image = None, sinogram_value = None, reconstruction_image = None):
    value = 255
    x1,y1 = start
    x2,y2 = end

    xi = 1
    yi = 1
    if x2 - x1 < 0:
        xi = -1
    if y2 - y1 < 0:
        yi = -1

    dx = np.abs(x2 - x1)
    dy = np.abs(y2 - y1)
    e = dx / 2
    sum = 0
    count = 0
    if dx == 0 or dy == 0:
        if reconstruction_image is None:
            s, c = draw_straight_line(image,oryg_image, start, end, value)
            sum += s
            count += c
        else:
            reconstruct_straight_line(reconstruction_image, sinogram_value, start, end)
    else:
        if reconstruction_image is None:
            image[y1, x1] = value
            sum += oryg_image[y1, x1]
            count += 1
        else:
            reconstruction_image[y1, x1] += sinogram_value
        if dx > dy:
            #OX priority
            for i in range(0, dx):
                x1 = x1 + xi
                e = e - dy
                if e < 0:
                    y1 = y1 + yi
                    e = e + dx
                if reconstruction_image is None:
                    image[y1, x1] = value
                    sum += oryg_image[y1, x1]
                    count += 1
                else:
                    reconstruction_image[y1, x1] += sinogram_value

        else:
            #OY priority
            for i in range(0, dy):
                y1 = y1 + yi
                e = e - dx
                if e < 0:
                    x1 = x1 + xi
                    e = e + dy
                if reconstruction_image is None:
                    image[y1, x1] = value
                    sum += oryg_image[y1, x1]
                    count += 1
                else:
                    reconstruction_image[y1, x1] += sinogram_value


    return (sum / (count+1) if reconstruction_image is None else None)

def reconstruct_straight_line(image, value, start, end):
    if start[0] == end[0]:
        image[:, start[0]] += value
    if end[1] == start[1]:
        image[start[1], :] += value

def draw_straight_line(image,oryg_image, start, end, value):
    #print("drawing straight line")
    if start[0] == end[0]:
        s = start[1] if start[1] < end[1] else end[1]
        e = end[1] if end[1] > start[1] else start[1]
        sum = np.sum(oryg_image[s:e, start[0]])
        count = e - s
        image[s:e, start[0]] = value

    if end[1] == start[1]:
        s = start[0] if start[0] < end[0] else end[0]
        e = end[0] if end[0] > start[0] else start[0]
        sum = np.sum(oryg_image[start[1], s:e])
        count = e - s

        image[start[1], s:e] = value
    return sum, count

def calculatePosition(angle, radius, circle_middle = (0,0)):
    rad = angle / 180 * np.pi
    return (circle_middle[0] + (int)(radius*np.cos(rad)), circle_middle[1] - (int)(radius*np.sin(rad)))

'''Calculates position x1, y1 and changes these values according to the available space'''
def calculatePositionSafe(angle, diameter, circle_middle = (0,0)):
    x1, y1 = calculatePosition(angle, diameter // 2 - 1, circle_middle)
    return check_borders((x1, y1), diameter)

def check_value(val, size):
    if val <= 0:
        val = 1
    if val >=size:
        val = size -1
    return val

def check_borders(point, size):
    return (check_value(point[0], size),check_value(point[1], size))

def reconstruct_line(sinogram_value, reconstruction_image, start_point, end_point):
    bresenham(sinogram_value=sinogram_value, reconstruction_image=reconstruction_image, start=start_point, end=end_point)

#draws line between start_point and end_point or draws a diameter using the angle
def draw_line(oryg_image, image, angle = None, diameter = None, center = None, start_point = None, end_point = None):
    if start_point != None and end_point != None:
        return bresenham(start_point, end_point, oryg_image=oryg_image, image=image)

    if diameter == None:
        diameter = len(image)

    radius = diameter // 2 - 1

    if center == None:
        center = (radius, radius)

    x1, y1 = calculatePositionSafe(angle, diameter, center)

    x2 = diameter - x1
    y2 = diameter - y1
    return bresenham((x1, y1), (x2, y2), oryg_image=oryg_image, image=image)
