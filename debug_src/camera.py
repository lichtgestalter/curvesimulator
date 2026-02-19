import cv2

ULIHOME = 0
ASPIT = 1

capture = cv2.VideoCapture(ASPIT)
while True:
    ret, frame = capture.read()
    if not ret:
        break
    print("Image taken")
    cv2.imshow("some window name", frame)
    cv2.waitKey(0)