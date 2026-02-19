import cv2

ULIHOME = 0
ASPIT = 1

camera = cv2.VideoCapture(ULIHOME)
while True:
    ret, frame = camera.read()
    print(ret)
    cv2.imshow("some window name", frame)
    cv2.waitKey(0)