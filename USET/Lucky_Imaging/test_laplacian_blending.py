import cv2
import numpy as np,sys

shape = np.array([512, 512])

A = np.ones(shape) * 80
B = np.ones(shape) * 255

# Colorize inside each side
A[200:300, 0:100] = 255
B[200:300, 400:] = 0

# Define large image frames acting as co-spatial canvas, hosting image A and B down below
CA = np.zeros(2*shape)
CB = np.zeros(2*shape)

CA[0:shape[0], 0:shape[1]] = A
CB[shape[0]/2:shape[0]/2+shape[0], shape[1]/2:shape[1]/2+shape[1]] = B

cv2.imwrite('CA.jpg', CA)
cv2.imwrite('CB.jpg', CB)

A = CA
B = CB

# generate Gaussian pyramid for A
G = A.copy()
gpA = [G]
for i in range(6):
    G = cv2.pyrDown(G)
    gpA.append(G)

# generate Gaussian pyramid for B
G = B.copy()
gpB = [G]
for i in range(6):
    G = cv2.pyrDown(G)
    gpB.append(G)

# generate Laplacian Pyramid for A
lpA = [gpA[5]]
for i in range(5, 0, -1):
    GE = cv2.pyrUp(gpA[i])
    L = cv2.subtract(gpA[i-1], GE)
    lpA.append(L)

# generate Laplacian Pyramid for B
lpB = [gpB[5]]
for i in range(5, 0, -1):
    GE = cv2.pyrUp(gpB[i])
    L = cv2.subtract(gpB[i-1], GE)
    lpB.append(L)

# Now add left and right halves of images in each level
LS = []
for la, lb in zip(lpA, lpB):
    rows, cols = la.shape
    ls = np.hstack((la[:, 0:cols/2], lb[:, cols/2:]))
    LS.append(ls)

# Now reconstruct
ls_ = LS[0]
for i in range(1,6):
    ls_ = cv2.pyrUp(ls_)
    ls_ = cv2.add(ls_, LS[i])

# image with direct connecting each half
real = np.hstack((A[:, :cols/2], B[:, cols/2:]))

cv2.imwrite('Pyramid_blending.jpg', ls_)
#cv2.imwrite('Direct_blending.jpg', real)
