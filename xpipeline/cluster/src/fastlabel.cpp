// replacement for bwlabel, the tweek beeing it starts with list of
// pixel coordinates instead of black-white logical maps, this should
// save computational time (at the price of complexity) because the
// former is a factor 100 smaller array than the latter
//
// inputs:
// coords - two column array with the x and y index of each pixel
// coordDim - two element vector with size of the x and y columns of the initial TF map
// outputs:
// labelList - vector with the same number of rows as coords
//             containing the label (cluster number) of each pixel in
//             the coords input

// author: michal.was@ligo.org

#include "fastlabel.h"

void fastlabel(const int nPixels, const double *coords, const double *coordDim, const int nNeighbors, double *labelList)
{
  int labelCounter = 0;
  queue<int> pixelQueue;

  vector<int> neighboorsX(nNeighbors);
  vector<int> neighboorsY(nNeighbors);

  if ( nNeighbors == 8 )
    {
      // 8 connectivity, a 3x3 square, without the middle pixel
      int iNeigh=0;
      for(int ix=-1; ix<=1 ; ix++)
        for(int iy=-1; iy<=1; iy++)
          {
            if ( ix==0 && iy==0 ) continue; //skip middle pixel
            neighboorsX[iNeigh] = ix;
            neighboorsY[iNeigh] = iy;
            iNeigh++;
          }
    }
  else if (nNeighbors == 24)
    {
      // 24 connectivity, a 5x5 square, without the middle pixel
      int iNeigh=0;
      for(int ix=-2; ix<=2 ; ix++)
        for(int iy=-2; iy<=2; iy++)
          {
            if ( ix==0 && iy==0 ) continue; //skip middle pixel
            neighboorsX[iNeigh] = ix;
            neighboorsY[iNeigh] = iy;
            iNeigh++;
          }
    }
  else if (nNeighbors == 48)
    {
      // 48 connectivity, a 7x7 square, without the middle pixel
      int iNeigh=0;
      for(int ix=-3; ix<=3 ; ix++)
        for(int iy=-3; iy<=3; iy++)
          {
            if ( ix==0 && iy==0 ) continue; //skip middle pixel
            neighboorsX[iNeigh] = ix;
            neighboorsY[iNeigh] = iy;
            iNeigh++;
          }
    }
  else if (nNeighbors == 80)
    {
      // 80 connectivity, a 9x9 square, without the middle pixel
      int iNeigh=0;
      for(int ix=-4; ix<=4 ; ix++)
        for(int iy=-4; iy<=4; iy++)
          {
            if ( ix==0 && iy==0 ) continue; //skip middle pixel
            neighboorsX[iNeigh] = ix;
            neighboorsY[iNeigh] = iy;
            iNeigh++;
          }
    }
  else
    {
      std::cout << "Error the connectivity number is not recognized, only 8, 24, 48 and 80 are currently supported " << std::endl;
      return;
    }


  // record conversion from pixel number to linear index in TF map
  vector<int> pixToLinIdx(nPixels);
  for(int iPix=0; iPix < nPixels; iPix++)
    {
      pixToLinIdx[iPix] = (coords[iPix] - 1)*coordDim[1] + coords[iPix + nPixels] - 1;
      // check that pixels are sorted (required for the algorithm to work)
      if(iPix>0 && pixToLinIdx[iPix-1] >= pixToLinIdx[iPix])
        std::cout << "Pixels are not sorted according to the linear index in the TF map"  << std::endl;
    }
  
  // create array of graph edges (links to neighboors)
  vector<int> edges(nPixels*nNeighbors);
  for(int iPix=0; iPix < nPixels; iPix++)
    {
      for( int iNeigh=0; iNeigh < nNeighbors; iNeigh++)
        {
          int newX = coords[iPix] - 1 + neighboorsX[iNeigh];
          int newY = coords[iPix + nPixels] - 1 + neighboorsY[iNeigh];
          // record linear index of neighboors or -1 if the neighboor
          // would be out of the TF map
          if ( newX >= 0 && newX < coordDim[0] && newY >= 0 && newY < coordDim[1])
            edges[iPix*nNeighbors + iNeigh] = newX*coordDim[1] + newY;
          else
            edges[iPix*nNeighbors + iNeigh] = -1;

        }
    }

  // replace linear indexes by pixel number in edge array and prune
  // edges that go to non black pixels
  for( int iNeigh=0; iNeigh < nNeighbors; iNeigh++)
    {
      int targetPix = 0;
      for(int iPix=0; iPix < nPixels; iPix++)
        {
          if (edges[iPix*nNeighbors + iNeigh] < 0)
            continue;
          while (edges[iPix*nNeighbors + iNeigh] > pixToLinIdx[targetPix] &&
                 targetPix < nPixels)
            targetPix++;
          if (targetPix < nPixels && 
              edges[iPix*nNeighbors + iNeigh] == pixToLinIdx[targetPix])
            edges[iPix*nNeighbors + iNeigh] = targetPix;
          else
            edges[iPix*nNeighbors + iNeigh] = -1;
        }
    }

  // loop over each pixel and start a new cluster (label) for each
  // pixel that is not labelled yet
  for(int iPix=0; iPix < nPixels; iPix++)
    {
      // skip pixel if already in cluster
      if( labelList[iPix] > 0 ) continue;

      // start new cluster (label)
      labelCounter++;
      labelList[iPix] = labelCounter;
      pixelQueue.push(iPix);
      // make a breath first search of the cluster
      while (!pixelQueue.empty())
        {
          int curPix = pixelQueue.front();
          pixelQueue.pop();

          for( int iNeigh=0; iNeigh < nNeighbors; iNeigh++)
            {
              if (edges[curPix*nNeighbors + iNeigh] >= 0)
                {
                  int newPix = edges[curPix*nNeighbors + iNeigh];
                  if( labelList[newPix] == 0)
                    {
                      labelList[newPix] = labelCounter;
                      pixelQueue.push(newPix);
                    }
                  else if ( labelList[newPix] != labelCounter )
                    std::cout << "Error : Trying to connect to an old cluster, this means the cluster construction algorithm is broken!" << std::endl;
                }

            }// loop over neighboors
        }// loop over pixels in cluster
    }// loop over all pixels

  return;
}
