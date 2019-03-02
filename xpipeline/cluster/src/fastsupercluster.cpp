// It is assumed that the input format is: [x y dx dy significance]
// For purposes of algorithm representation we interpret
// significance as height, so superclustering is just looking for
// bounding boxes that are not covered by any other
// For two overlapping clusters with equal significance the cluster
// which comes earlier in the list is kept

// same as fastquadraticsupercluster.cpp but with better complexity,
// should scale for large amounts of events (i.e. superclustering over
// hundreds of sky positions)

// Authors:
//
//    Dorota Was
//    Michal Was    michal.was@ligo.org

#include "fastsupercluster.h"

class Point{
public: 
  double YCor;
  int idx;
  
  Point(const double y, const int i) {
    YCor = y;
    idx = i;
  }
  Point() {
    YCor = 0;
    idx = -1; 
  }
  
  bool operator < (const Point & A) const {
    if(YCor != A.YCor)
      return YCor < A.YCor;
    return idx < A.idx;
  }
};

bool intersectSeg(double begA, double deltaA, double begB, double deltaB) { // Segments are parallel
  if((begA <= begB) && begB < (begA + deltaA))
    return true;
  if((begB <= begA) && begA < (begB + deltaB)) 
    return true; 
  return false;
}

bool intersectRec(const int N, const double* matrice, int i, int j) {
  bool result = true;
  if(!intersectSeg(matrice[i], matrice[i+2*N], matrice[j], matrice[j+2*N])) // Checking x-coordinate
    result = false;
  if(!intersectSeg(matrice[i+N], matrice[i+3*N], matrice[j+N], matrice[j+3*N])) // Checking y-coordinate
    result = false;
  return result;
}

void fastsupercluster(const int N, const double* rectangles, double* uncoveredMask) { 
// returns an array T, such that T[i] = true, if the i-th rectangle is uncovered, = false otherwise 
  // do: delete[] T;, after use
   set<double> XCor;
  XCor.clear();
  for(int i=0; i<N; i++) {
    uncoveredMask[i] = true;
    XCor.insert(rectangles[i]);
    XCor.insert(rectangles[i] + rectangles[i+2*N]);
  }
  int size = XCor.size();
  double* XCorT = new double[size];
  int idx = 0;
  for(set<double>::iterator i = XCor.begin(); i!=XCor.end(); i++) {
    XCorT[idx++] = *i;
  }
  list<int>* Indexes = new list<int>[size];
  for(int i=0; i<N; i++) {
    double* p = lower_bound(XCorT, XCorT + size, rectangles[i]);
    idx = p - XCorT;
    Indexes[idx].push_back(i); // so that inserting new segments will be done after removing
    p = lower_bound(XCorT, XCorT + size, rectangles[i] + rectangles[i + 2*N]);
    idx = p - XCorT;
    Indexes[idx].push_front(i); // so that removing segments will be done after inserting the new ones 
  }
  
  set<Point> begins;
  for(int i=0; i<size; i++) {
    for(list<int>::iterator iter = Indexes[i].begin(); iter!=Indexes[i].end(); iter++) {
      if(rectangles[*iter] + rectangles[*iter + 2*N] <= XCorT[i]) { // removes the segment (rectangles[*iter+N] , rectangles[*iter+N] + rectangles[*iter+3*N])
	set<Point>::iterator temp = begins.find(Point(rectangles[*iter+N], *iter));
	if(temp == begins.end())
	  printf("Error. Cannot remove a rectangle. Numerical problem?\n");
	else
	  begins.erase(temp);
      }
      else { // adds the segment (rectangles[*iter+N] , rectangles[*iter+N] + rectangles[*iter+3*N])
	set<Point>::iterator begIt = begins.begin();
	set<Point>::iterator endIt = begins.upper_bound(Point(rectangles[*iter+N] + rectangles[*iter+3*N], -1));
	
	for(set<Point>::iterator pointsIt = begIt; pointsIt != endIt; pointsIt++) {
	  if(intersectSeg(rectangles[*iter + N], rectangles[*iter + 3*N], rectangles[pointsIt->idx + N], rectangles[pointsIt->idx + 3*N])) {
	    if(rectangles[*iter + 4*N] > rectangles[pointsIt->idx + 4*N])
	      uncoveredMask[pointsIt->idx] = false;
	    else if(rectangles[*iter + 4*N] < rectangles[pointsIt->idx + 4*N])
	      uncoveredMask[*iter] = false;
	    else if(*iter < pointsIt->idx)
	      uncoveredMask[pointsIt->idx] = false;
	    else
	      uncoveredMask[*iter] = false;
	  }
	}
	begins.insert(Point(rectangles[*iter+N], *iter));
      }
    }
  }
  delete[] Indexes;
  delete[] XCorT;
}
