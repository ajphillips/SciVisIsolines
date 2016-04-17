/*=========================================================================

    This project uses the Visualization Toolkit library and helper functions
    from TriangleList.h to render an isosurface. Implementation by Andrew Phillips.

=========================================================================*/

#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

int main()
{
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("Isolines.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

	int vert[4][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 } };
	int ptIdx[4] = { 0, 0, 0, 0 };
	float endPoints[4][2];

    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

	//Algorithm for drawing lines between endpoints on cells
	for (int y = 0; y < 49; y++){
		for (int x = 0; x < 49; x++){
			int pts = 0;

			endPoints[0][0] = 0;
			endPoints[0][1] = 0;
			endPoints[1][0] = 0;
			endPoints[1][1] = 0;
			endPoints[2][0] = 0;
			endPoints[2][1] = 0;
			endPoints[3][0] = 0;
			endPoints[3][1] = 0;

			vert[0][0] = x;
			vert[0][1] = y;
			vert[1][0] = x + 1;
			vert[1][1] = y;
			vert[2][0] = x;
			vert[2][1] = y + 1;
			vert[3][0] = x + 1;
			vert[3][1] = y + 1;

			ptIdx[0] = GetPointIndex(vert[0], dims);
			ptIdx[1] = GetPointIndex(vert[1], dims);
			ptIdx[2] = GetPointIndex(vert[2], dims);
			ptIdx[3] = GetPointIndex(vert[3], dims);

			if (F[ptIdx[0]] < 3.2000){
				pts += 1;
			}
			if (F[ptIdx[1]] < 3.2000){
				pts += 10;
			}
			if (F[ptIdx[2]] < 3.2000){
				pts += 100;
			}
			if (F[ptIdx[3]] < 3.2000){
				pts += 1000;
			}

			// ******************************************************************************************
			/*
			I split my cases up into these statements instead of using a look up table.

			Because the function to draw lines just needs the endpoints of the line, I can combine cases
			and do not have to worry about measurements like color or value
			I was going to convert the following if statements into a switch but I felt it looked even
			more cluttered with that so I left it as is.
			*/
			// ******************************************************************************************

			if (pts == 0){
				//Case 0
				//Do nothing
			}
			else if (pts == 1 || pts == 1110){
				//Cases 1 and 14
				endPoints[0][0] = X[x] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[1]] - F[ptIdx[0]]))*(X[x + 1] - X[x]);
				endPoints[0][1] = Y[y];
				endPoints[1][0] = X[x];
				endPoints[1][1] = Y[y] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[2]] - F[ptIdx[0]]))*(Y[y + 1] - Y[y]);

				sl.AddSegment(endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]);
			}
			else if (pts == 10 || pts == 1101){
				//Cases 2 and 13
				endPoints[0][0] = X[x] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[1]] - F[ptIdx[0]]))*(X[x + 1] - X[x]);
				endPoints[0][1] = Y[y];
				endPoints[1][0] = X[x + 1];
				endPoints[1][1] = Y[y] + ((3.2 - F[ptIdx[1]]) / (F[ptIdx[3]] - F[ptIdx[1]]))*(Y[y + 1] - Y[y]);

				sl.AddSegment(endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]);
			}
			else if (pts == 11 || pts == 1100){
				//Cases 3 and 12
				endPoints[0][0] = X[x];
				endPoints[0][1] = Y[y] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[2]] - F[ptIdx[0]]))*(Y[y + 1] - Y[y]);
				endPoints[1][0] = X[x + 1];
				endPoints[1][1] = Y[y] + ((3.2 - F[ptIdx[1]]) / (F[ptIdx[3]] - F[ptIdx[1]]))*(Y[y + 1] - Y[y]);

				sl.AddSegment(endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]);
			}
			else if (pts == 100 || pts == 1011){
				//Cases 4 and 11
				endPoints[0][0] = X[x];
				endPoints[0][1] = Y[y] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[2]] - F[ptIdx[0]]))*(Y[y + 1] - Y[y]);
				endPoints[1][0] = X[x] + ((3.2 - F[ptIdx[2]]) / (F[ptIdx[3]] - F[ptIdx[2]]))*(X[x + 1] - X[x]);
				endPoints[1][1] = Y[y + 1];

				sl.AddSegment(endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]);
			}
			else if (pts == 101 || pts == 1010){
				//Cases 5 and 10
				endPoints[0][0] = X[x] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[1]] - F[ptIdx[0]]))*(X[x + 1] - X[x]);
				endPoints[0][1] = Y[y];
				endPoints[1][0] = X[x] + ((3.2 - F[ptIdx[2]]) / (F[ptIdx[3]] - F[ptIdx[2]]))*(X[x + 1] - X[x]);
				endPoints[1][1] = Y[y + 1];

				sl.AddSegment(endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]);
			}
			else if (pts == 110){
				//Case 6
				endPoints[0][0] = X[x] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[1]] - F[ptIdx[0]]))*(X[x + 1] - X[x]);
				endPoints[0][1] = Y[y];
				endPoints[1][0] = X[x + 1];
				endPoints[1][1] = Y[y] + ((3.2 - F[ptIdx[1]]) / (F[ptIdx[3]] - F[ptIdx[1]]))*(Y[y + 1] - Y[y]);
				endPoints[2][0] = X[x] + ((3.2 - F[ptIdx[2]]) / (F[ptIdx[3]] - F[ptIdx[2]]))*(X[x + 1] - X[x]);
				endPoints[2][1] = Y[y + 1];
				endPoints[3][0] = X[x];
				endPoints[3][1] = Y[y] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[2]] - F[ptIdx[0]]))*(Y[y + 1] - Y[y]);

				sl.AddSegment(endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]);
				sl.AddSegment(endPoints[2][0], endPoints[2][1], endPoints[3][0], endPoints[3][1]);
			}
			else if (pts == 111 || pts == 1000){
				//Cases 7 and 8
				endPoints[0][0] = X[x] + ((3.2 - F[ptIdx[2]]) / (F[ptIdx[3]] - F[ptIdx[2]]))*(X[x + 1] - X[x]);
				endPoints[0][1] = Y[y + 1];
				endPoints[1][0] = X[x + 1];
				endPoints[1][1] = Y[y] + ((3.2 - F[ptIdx[1]]) / (F[ptIdx[3]] - F[ptIdx[1]]))*(Y[y + 1] - Y[y]);

				sl.AddSegment(endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]);
			}
			else if (pts == 1001){
				//Case 9
				endPoints[0][0] = X[x] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[1]] - F[ptIdx[0]]))*(X[x + 1] - X[x]);
				endPoints[0][1] = Y[y];
				endPoints[1][0] = X[x];
				endPoints[1][1] = Y[y] + ((3.2 - F[ptIdx[0]]) / (F[ptIdx[2]] - F[ptIdx[0]]))*(Y[y + 1] - Y[y]);
				endPoints[2][0] = X[x] + ((3.2 - F[ptIdx[2]]) / (F[ptIdx[3]] - F[ptIdx[2]]))*(X[x + 1] - X[x]);
				endPoints[2][1] = Y[y + 1];
				endPoints[3][0] = X[x + 1];
				endPoints[3][1] = Y[y] + ((3.2 - F[ptIdx[1]]) / (F[ptIdx[3]] - F[ptIdx[1]]))*(Y[y + 1] - Y[y]);

				sl.AddSegment(endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]);
				sl.AddSegment(endPoints[2][0], endPoints[2][1], endPoints[3][0], endPoints[3][1]);
			}
			else if (pts == 1111){
				//Case 15
				//Do nothing
			}
			else{
				//Do nothing
			}
			//End of inner loop
		}
		//End of outer loop
	}
	//End of algorithm

    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
