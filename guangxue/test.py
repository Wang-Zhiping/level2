# import cv2
# import numpy as np

# if __name__ == '__main__':
#     img = cv2.imread("D:\Code\PY\guangxue\p.png") 
#     # cv2.imshow("Image", img) 
#     # cv2.waitKey(0)

#     lower = np.array([1,120,30], np.uint8)
#     upper = np.array([255,200,100], np.uint8)
	
#     height, width = img.shape[:2]
    
#     abc = cv2.cvtColor(img, cv2.COLOR_BGR2YCR_CB)
#     cv2.imshow("Image", cv2.resize(abc, (int(width/2), int(height/2)), interpolation=cv2.INTER_CUBIC)) 
#     cv2.waitKey(0)
#     mask = cv2.inRange(abc, lower, upper)
#     cv2.imshow("Image", cv2.resize(mask, (int(width/2), int(height/2)), interpolation=cv2.INTER_CUBIC)) 
#     cv2.waitKey(0)

#     kernel = np.ones((57,57),np.uint8)
#     opening = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
#     skin = cv2.bitwise_and(img, img, mask=opening)
#     cv2.imwrite("D:\Code\PY\guangxue\{}skin.jpg".format(1),skin)

    # cv2.imshow("Image", cv2.resize(skin, (int(width/2), int(height/2)), interpolation=cv2.INTER_CUBIC)) 
    # cv2.waitKey(0)


# import cv2
# import numpy as np

# if __name__ == '__main__':
#     img = cv2.imread("D:\Code\PY\guangxue\p.png") 
#     # cv2.imshow("Image", img) 
#     # cv2.waitKey(0)

#     a,b,c,d,e,f=80,0,0,255,0,0
#     z = range(0,255,50)

#     for b in z:
#       for c in z:
#         for e in range(b,255,50):
#          for f in range(c,255,50):
#             lower = np.array([a,b,c], np.uint8)
#             upper = np.array([d,e,f], np.uint8)
            
#             abc = cv2.cvtColor(img, cv2.COLOR_BGR2YCR_CB)
#             # cv2.imshow("Image", abc) 
#             # cv2.waitKey(0)
#             mask = cv2.inRange(abc, lower, upper)

#             kernel = np.ones((57,57),np.uint8)
#             opening = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
#             skin = cv2.bitwise_and(img, img, mask=opening)
#             cv2.imwrite("D:\Code\PY\guangxue\{},{},{}_{},{},{}skin.jpg".format(a,b,c,d,e,f),skin)
#             # # cv2.imshow("Image", skin) 
#             # cv2.waitKey(0)

import cv2
import numpy as np
import matplotlib.pyplot as plt
def gray(i,img):
    img_gray = cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)
    cv2.imwrite("gray.jpg",img_gray)
    cv2.imshow("img_gray{}".format(i),cv2.resize(img_gray, (int(width/2), int(height*10)), interpolation=cv2.INTER_CUBIC))



if __name__ == '__main__':
    img = cv2.imread("D:\Code\PY\guangxue\pp.png") 
    height, width = img.shape[:2]

    # cv2.imshow("Image",cv2.resize(img, (int(width/2), int(height*10)), interpolation=cv2.INTER_CUBIC))
    gray(1,img)
    cv2.waitKey(0)

#    print(cv2.cvtColor(img,cv2.COLOR_RGB2GRAY))

    #亮度分布
    plt.figure(num=1,figsize=(7,5))
    x = np.linspace(0,len(img[0]),len(img[0]))
    l=np.zeros(len(img[0]))

    for i in range(len(img)):
        plt.plot(cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)[i],alpha=0.05)
        for j in range(len(img[0])):
            l[j] += cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)[i][j]
    l=l/len(img)

    N=10
    for i in range(len(l)):
        l[i] = np.mean(l[i:i+N])

    plt.plot(x,l,color="Tab:blue",alpha=0.7)

    ml = []
    mx = []
    for i in range(1,len(l)-1):
        if (l[i-1]<=l[i])and(l[i+1]<=l[i]):
            ml.append(l[i])
            mx.append(x[i])
    plt.scatter(mx,ml,alpha=0.2,color="Tab:orange",edgecolors="Tab:orange")







    #亮度梯度
    ll=l.copy()
    for i in range(len(l)-1):
        ll[i+1]=l[i+1]-l[i]
    ll[0] = 0

    lll=ll.copy()
    for i in range(len(l)-1):
        lll[i+1]=ll[i+1]-ll[i]
    lll[0] = 0


    plt.figure(num=2,figsize=(15,5))
    plt.ylim(-7.5,7.5)
    plt.plot(x,ll,linewidth=0.5)
    plt.xlabel("dy")
    plt.plot([0,len(img[0])],[0,0],linewidth=0.5)

    mll = []
    mxx = []
    mml = []
    for i in range(1,len(ll)-1):
        if (abs(ll[i-1])>=abs(ll[i]))and(abs(ll[i+1])>=abs(ll[i])):
            mml.append(l[i])
            mll.append(ll[i])
            mxx.append(x[i])
    plt.scatter(mxx,mll,alpha=0.2,color="Tab:orange",edgecolors="Tab:orange")
    
    #交集
    _x = np.intersect1d(mx,mxx)
    _x = np.delete(_x,np.where(_x<710))
    _x = np.delete(_x,np.where(_x>910))
    x1 = np.delete(np.array(mx),np.where(np.array(mx) > 710))
    x2 = np.delete(np.array(mx),np.where(np.array(mx) < 970))
    _x = np.union1d(x1,_x)
    _x = np.union1d(_x,x2)
    
    for i in range(100):
        if (_x[i+1]-_x[i])<2:
            _x = np.delete(_x,np.where(_x==_x[i]))
    order=[]
    for i in range(0,len(x)):
        sum = 0
        for j in range(0,len(_x)):
            if _x[j]==x[i]:
                sum += 1
        if sum==0:
            order.append(i)
    print(len(x)-len(order)-len(_x))
    _l = np.delete(l,order)
    plt.figure(num=3,figsize=(15,5))
    plt.plot(x,l,linewidth=2,alpha=0.5)
    plt.scatter(_x,_l,alpha=0.4,color="Tab:orange",edgecolors="Tab:orange")


    a = np.where(_l<100)
    _l = np.delete(_l,a)
    _x = np.delete(_x,a)
    plt.figure(num=4,figsize=(15,5))
    plt.plot(x,l,linewidth=2,alpha=0.5)
    plt.scatter(_x,_l,alpha=0.4,color="Tab:orange",edgecolors="Tab:orange")

    # plt.scatter(mxx,mml,alpha=0.4,color="Tab:green",edgecolors="Tab:green")


    # plt.figure(num=3,figsize=(15,5))
    # plt.ylim(-100,100)
    # plt.plot(x,lll,linewidth=0.5)
    # plt.plot([0,len(img[0])],[0,0],linewidth=0.5)
    # plt.xlabel("ddy")
    
    plt.show() 

    plt.plot(np.linspace(1,len(_x),len(_x)),_x)
    plt.show()

    print(_x)