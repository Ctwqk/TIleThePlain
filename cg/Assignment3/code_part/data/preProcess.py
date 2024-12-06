import cv2
import numpy as np
jpg_ = cv2.imread("3333.jpg")

jpg = cv2.cvtColor(jpg_,cv2.COLOR_BGR2GRAY)

cv2.imshow("jpg_",jpg_)
#optm_thshd, res = cv2.threshold(jpg,0,1,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
#cv2.imshow("res",res)

#mask = res.astype(bool)

#alpha = np.ones(jpg.shape,dtype = jpg.dtype)*255

#jpg_ = cv2.merge((jpg_,alpha))
#jpg_[mask,3] = 0

#cv2.imshow("final",jpg_)

#cv2.imwrite("tile.png", jpg_)


cv2.waitKey()
cv2.destroyAllWindows()

#processed = cv2.imread("tile.png", cv2.IMREAD_UNCHANGED)