import cv2
camera = cv2.VideoCapture(0)
while True:
    ret, frame = camera.read()
    print(ret)
    cv2.imshow("some window name", frame)
    cv2.waitKey(0)