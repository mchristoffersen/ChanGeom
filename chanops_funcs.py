# Michael Christoffersen
# python 2.7

# File Summary
# This file contains all of the functions that were written to be used with the chanops.py script
# The functions are as follows:
# A binary hit-or-miss function, an endpoint detection function for a skeletonized binary image,
# a function that finds all of the pixels set to one in the 8 pixels surrounding each pixel in a
# list of input pixels, a function that takes two lists of pixel coordinates and checks if there
# are any neighboring pixels contained in the two lists, a function to determine if two points are
# diagonal to each other or not, and a function to trace up a single pixel channel skeleton while
# keeping track of distance.
#
# Note, all pixel lists in these scripts are in the format ([a,b,c],[d,e,f]) where (a,d) is a point,
# (b,e) is a point and so on. So to access the first coordinate pair in a list of points called "points"
# the indicies to access would be points[0][0] and points[1][0]

import numpy as np

# A binary pattern detection function, takes a numpy array representation of the image and
# a 3x3 numpy array of the desired pattern as inputs. The function will match
# the ones and zeroes in the input array and ignore all other numbers. Returns an array of the
# same dimensions as the input image array with ones where the pattern was detected and zeros
# everywhere else
def binhitmiss(binar,ptn,idx_check):
	(bx,by) = binar.shape
	idx_true = np.zeros((bx,by))
	idx1 = list(np.where(ptn==1))
	idx1[0] = idx1[0] - 1
	idx1[1] = idx1[1] - 1

	idx0 = list(np.where(ptn==0))
	idx0[0] = idx0[0] - 1
	idx0[1] = idx0[1] - 1

	for idx in range(len(idx_check[0])):
		i = idx_check[0][idx]
		j = idx_check[1][idx]
		onesmatch = True
		zerosmatch = True

		for idxx in range(len(idx1[0])):
			m = i + idx1[0][idxx]
			n = j + idx1[1][idxx]
			if m < 0 or n < 0 or m >= bx or n >= by:
				onesmatch = False
			elif binar[m][n] != 1:
				onesmatch = False

		for idxx in range(len(idx0[0])):
			m = i + idx0[0][idxx]
			n = j + idx0[1][idxx]
			if m < 0 or n < 0 or m >= bx or n >= by:
				zerosmatch = False
			elif binar[m][n] != 0:
				zerosmatch = False

		if onesmatch == True and zerosmatch == True:
			idx_true[i][j] = 1

	return idx_true

#shifted binary hit or miss
def binhitmiss_shift(binar,ptn,idx_check,shift):
	#shift = 1 is to the right, -1 to the left
	(bx,by) = binar.shape
	idx_true = np.zeros((bx,by))
	idx1 = list(np.where(ptn==1))
	idx1[0] = idx1[0] - 1
	idx1[1] = idx1[1] - 1

	idx0 = list(np.where(ptn==0))
	idx0[0] = idx0[0] - 1
	idx0[1] = idx0[1] - 1

	for idx in range(len(idx_check[0])):
		i = idx_check[0][idx]
		j = idx_check[1][idx]
		onesmatch = True
		zerosmatch = True

		for idxx in range(len(idx1[0])):
			m = i + idx1[0][idxx]
			n = j + idx1[1][idxx] + shift
			if m < 0 or n < 0 or m >= bx or n >= by:
				onesmatch = False
			elif binar[m][n] != 1:
				onesmatch = False

		for idxx in range(len(idx0[0])):
			m = i + idx0[0][idxx]
			n = j + idx0[1][idxx] + shift
			if m < 0 or n < 0 or m >= bx or n >= by:
				zerosmatch = False
			elif binar[m][n] != 0:
				zerosmatch = False

		if onesmatch == True and zerosmatch == True:
			idx_true[i][j] = 1

	return idx_true

# Function to find all of the endpoints in a channel skeleton. Takes in a binary numpy array
# where ones represent the channel skeleton and finds all the endpoints on the skeleton.
# Returns a list of the indicies of the endpoints.
def find_endpoints(imgar,idx_check):
	interv1 = np.array([[0,0,0],[0,1,0],[-1,-1,0]])
	interv2 = np.array([[0,0,0],[0,1,0],[0,-1,-1]])
	listends = [],[]
	for k in range(4):
		end1 = binhitmiss(imgar,np.rot90(interv1,k),idx_check)
		end2 = binhitmiss(imgar,np.rot90(interv2,k),idx_check)
		loc1 = list(np.where(end1==1))
		loc2 = list(np.where(end2==1))
		for item in range(len(loc1[0])):
			listends[0].append(loc1[0][item])
			listends[1].append(loc1[1][item])

		for idx in range(len(loc2[0])):
			listends[0].append(loc2[0][idx])
			listends[1].append(loc2[1][idx])

	return listends

# Finds the indicies of ones around input pixels
# Takes input "img", a binary image, and "locs", the list of points around
# which to search for ones.
def neighborones(img,locs):
	idxones = [],[]
	for idx in range(len(locs[0])):
		i = locs[0][idx]
		j = locs[1][idx]
		possib = [i-1, i-1, i-1, i, i, i+1, i+1, i+1],[j-1, j, j+1, j-1, j+1, j-1, j, j+1]
		for idxx in range(8):
			m = possib[0][idxx]
			n = possib[1][idxx]
			if img[m][n] == 1:
				idxones[0].append(m)
				idxones[1].append(n)

	return idxones

def thin(img,rcr):
	# either finite times or until image stops changing
	if(rcr == 0):
		while True:
			# thin img
			newimg = img[:]
			img = iter(img)
			if(np.array_equal(img,newimg)):
				break

	elif(rcr > 0):
		# code for finite thinning
		while (rcr > 0):
			# thin img
			img = iter(img)
			rcr = rcr - 1

	return img

def iter(img):
	#nhd variable is a 9 entry vector with pixel p in location 0 and starting with the east neighbor
	#put in in counterclockwise order:
	# 4 3 2
	# 5 p 1
	# 6 7 8
	#run g1+g2+g3
	#run g1+g2+g4
	#delete pixels
	locs = np.where(img == 1)
	newimg = np.copy(img)
	for i in range(len(locs[0])):
		nhd = getnhd(img,[locs[0][i], locs[1][i]])
		if(g1(nhd) and g2(nhd) and g3(nhd)):
			newimg[locs[0][i],locs[1][i]] = 0
	img = np.copy(newimg)
	for i in range(len(locs[0])):
		nhd = getnhd(img,[locs[0][i], locs[1][i]])
		if(g1(nhd) and g2(nhd) and g4(nhd)):
			newimg[locs[0][i],locs[1][i]] = 0

#	for i in range(len(locs[0])):
#		nhd = getnhd(img,[locs[0][i],locs[1][i]])
#		if(np.sum(nhd) <= 4):
#			continue
#		if(not (nhd[1] and nhd[3] and nhd[5] and nhd[7])):
#			newimg[locs[0][i],locs[1][i]] = 0

	return newimg

def g1(nhd):
	#g1 code
	xh = 0
	for i in [1,2,3,4]:
		xh = xh + int(not nhd[2*i-1] and (nhd[2*i] or nhd[(2*i+1)%8]))
	return (xh == 1)

def g2(nhd):
	#g2 code
	n1 = 0
	n2 = 0
	for i in [1,2,3,4]:
		n1 = n1 + int(nhd[2*i-1] or nhd[2*i])
	for i in [1,2,3,4]:
		n2 = n2 + int(nhd[2*i] or nhd[(2*i+1)%8])
	return (2 <= min(n1,n2) and min(n1,n2) <= 3)

def g3(nhd):
	#g3 code
	return (not ((nhd[2] or nhd[3] or (not nhd[8])) and nhd[1]))

def g4(nhd):
	#g4 code
	return (not ((nhd[6] or nhd[7] or (not nhd[4])) and nhd[5]))

def getnhd(img,loc):
	# img is a binary array
	# loc is a tuple or list in the format (row,col)
	lr = loc[0]
	lc = loc[1]
	unrl = np.zeros(9)
	unrl[0] = img[lr,lc]

	try:
		unrl[1] = img[lr,lc+1]
	except IndexError:
		unrl[1] = 0

	try:
		unrl[2] = img[lr-1,lc+1]
		if(lr == 0):
			unrl[2] = 0
	except IndexError:
		unrl[2] = 0

	try:
		unrl[3] = img[lr-1,lc]
		if(lr == 0):
			unrl[3] = 0
	except IndexError:
		unrl[3] = 0

	try:
		unrl[4] = img[lr-1,lc-1]
		if(lr == 0 or lc == 0):
			unrl[4] = 0
	except IndexError:
		unrl[4] = 0

	try:
		unrl[5] = img[lr,lc-1]
		if(lc == 0):
			unrl[5] = 0
	except IndexError:
		unrl[5] = 0

	try:
		unrl[6] = img[lr+1,lc-1]
		if(lc == 0):
			unrl[6] = 0
	except IndexError:
		unrl[6] = 0

	try:
		unrl[7] = img[lr+1,lc]
	except IndexError:
		unrl[7] = 0

	try:
		unrl[8] = img[lr+1,lc+1]
	except IndexError:
		unrl[8] = 0
	#print(unrl)
	return unrl



# Takes an input of two lists of points and compares the two lists, returning the points
# in the list "check" that are neighbors to points in the list "fixed"
def findneighbors(fixed,check):
	are_nabes = False
	nabe_locs = [],[]
	lenfix = len(fixed[0])
	lenchk = len(check[0])
	idx = 0
	while idx < lenchk:
		i = check[0][idx]
		j = check[1][idx]
		for idxx in range(lenfix):
			m = fixed[0][idxx]
			n = fixed[1][idxx]
			if abs(m-i) < 2 and abs(n-j) < 2:
				nabe_locs[0].append(i)
				nabe_locs[1].append(j)
				are_nabes = True
				del check[0][idx]
				del check[1][idx]
				lenchk = len(check[0])
		idx = idx + 1

	return nabe_locs,are_nabes,check

# Function determines if two points are diagonal or not. Takes two one point long lists
# of points and returns True if the two points
# are diagonal and False if they are on the same x or y axis.
def isdiag(pt1,pt2):
	diag = False
	i = pt1[0][0]
	j = pt1[1][0]
	m = pt2[0][0]
	n = pt2[1][0]
	ydif = i - m
	xdif = j - n
	if xdif != 0 and ydif != 0:
		diag = True
	return diag

# Function written specifically to iteratively trace up a river skeleton
# It takes in the current point on the skeleton, the point behind it in the skeleton trace,
# a list of the points in the skeleton, and then finds
# the next point in the skeleton and calculates a new total length of the skeleton.
# See the chanops.py script for proper usage.
def rivcount(curpt,laspt,img,distot):
	nabe_idxs = neighborones(img,curpt)
	nexpt = [0],[0]

	if laspt[0][0]!=nabe_idxs[0][0] or laspt[1][0]!=nabe_idxs[1][0]:
		nexpt[0][0] = nabe_idxs[0][0]
		nexpt[1][0] = nabe_idxs[1][0]
	elif laspt[0][0]!=nabe_idxs[0][1] or laspt[1][0]!=nabe_idxs[1][1]:
		nexpt[0][0] = nabe_idxs[0][1]
		nexpt[1][0] = nabe_idxs[1][1]
	else:
		print("Point Tracing Error")

	if isdiag(curpt,nexpt):
		distot = distot + 2**(0.5)
	else:
		distot = distot + 1

	return [nexpt,distot]
