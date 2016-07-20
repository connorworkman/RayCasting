/*============================================================
Author: Connor Workman

The idea here is fairly simple.  We project rays from a 
viewing point (a "camera") through a 3D data set.  
We take samplings of the values stored at several points
in the 3D data set and use these values to
assign a color to each ray (each pixel).  The result
is a 2D image that captures the size and color of the 3D data
from a particular viewpoint.

Notes:
* Has a curved "frustum" (near and far distance used for each ray),
but doesn't affect this dataset
* Uses brute force intersection evaluation instead of intelligently
determining cell intersections. (very slow at high sample rate/resolution)

============================================================*/

#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <float.h>
#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#define PIXEL_WIDTH 500
#define PIXEL_HEIGHT 500
#define NUM_SAMPLES 500
#define NUM_PICS 35


struct Camera
{
    double nearPlane, farPlane;
    double viewAngle;
    double viewPosition[3];
    double viewFocus[3];
    double viewUp[3];
};

struct ColorTransfer
{
    double minimum;
    double maximum;
    int numBins;
    unsigned char *colors;  // size is 3*numBins
    double *opacities; // size is numBins

    void ApplyColorTransfer(double value, unsigned char *RGB, double &opacity)
    {
        int bin;
        double bin_size = (maximum-minimum)/numBins;
        for (int i = 0; i < numBins; i++)
        {
            if (value > maximum)
            {
                bin = 255;
                RGB[0] = colors[3*bin+0];
                RGB[1] = colors[3*bin+1];
                RGB[2] = colors[3*bin+2];
                opacity = opacities[bin];
                return;
            }
            else if (value <= (minimum + bin_size*(i+1)))
            {
                bin = i;        
                RGB[0] = colors[3*bin+0];
                RGB[1] = colors[3*bin+1];
                RGB[2] = colors[3*bin+2];
                opacity = opacities[bin];
                return;
            }
        } 
    }
};

ColorTransfer
SetupColorTransfer(void)
{
    int  i;
    ColorTransfer rv;
    rv.minimum = 10;
    rv.maximum = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char opacity[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; //dictionary for mapping character opacity
    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = opacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = {71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76};
    double controlPointPositions[numControlPoints] = {0.0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0};
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
        if (end >= rv.numBins)
            end = rv.numBins-1;

        //set rv.colors
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.0)
                continue;
            for (int k = 0 ; k < 3 ; k++)
            {
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k]) + controlPointColors[3*i+k];
            }
        }
    }
    return rv;
}

/* SetupCamera returns a camera struct with attributes filled in. */
/* Nothing interesting here. */
Camera
SetupCamera(void)
{
    Camera rv;
    rv.viewFocus[0] = 0;
    rv.viewFocus[1] = 0;
    rv.viewFocus[2] = 0;
    rv.viewUp[0] = 0;
    rv.viewUp[1] = 1;
    rv.viewUp[2] = 0;
    rv.viewAngle = 30;
    rv.nearPlane = 7.5e+7;
    rv.farPlane = 1.4e+8;
    rv.viewPosition[0] = -8.25e+7;
    rv.viewPosition[1] = -1.45e+7;
    rv.viewPosition[2] = 3.35e+7;
    return rv;
}

int
GetPointIndex(const int x, const int y, const int z, const int *dims)
{
    return z*dims[0]*dims[1] + y*dims[0] + x;
}

int 
GetPointIndex(const int *idx, const int *dims)
{
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
}

/* This function takes the main components of the data set and produces a scalar value for a single point in the data */
/* In terms of linear algebra, this still isn't the interesting part. */
double EvaluateFieldAtLocation(const int *dims, const float *X, const float *Y, const float *Z, const float *F, const double *pt)
{
    int x, y, z;
    x = y = z = -1;
    for(int i = 0; i < dims[0]; i++)
    {
        if(pt[0] <= X[i])
        {
            x = i-1;
            for(int j=0; j < dims[1]; j++)
            {
                if(pt[1] <= Y[j])
                {
                    y = j-1;
                    for(int k = 0; k < dims[2]; k++)
                    {
                        if(pt[2] <= Z[k])
                        {
                            z = k-1;
                            break;
                        }
                    }
                    break;
                }
            }
            break;
        }
    }
    
    if (x == -1 || y == -1 || z == -1)
        return 0.0;

    /* Here we do some linear interpolation  to find approximate point value*/
    double x_prop = (pt[0] - X[x]) / (X[x+1] - X[x]);
    double x_0 = F[GetPointIndex(x,y,z, dims)]*(1-x_prop)
        +F[GetPointIndex(x+1,y,z, dims)]*x_prop;
    double x_1 = F[GetPointIndex(x,y,z+1, dims)]*(1-x_prop)
        +F[GetPointIndex(x+1,y,z+1, dims)]*x_prop;
    double x_2 = F[GetPointIndex(x,y+1,z, dims)]*(1-x_prop)
        +F[GetPointIndex(x+1,y+1,z, dims)]*x_prop;
    double x_3 = F[GetPointIndex(x,y+1,z+1, dims)]*(1-x_prop)
        +F[GetPointIndex(x+1,y+1,z+1, dims)]*x_prop;
    
    double y_prop = (pt[1] - Y[y]) / (Y[y+1] - Y[y]);
    double y_0 = x_0*(1-y_prop)+x_2*y_prop;
    double y_1 = x_1*(1-y_prop)+x_3*y_prop;

    double z_prop = (pt[2] - Z[z]) / (Z[z+1] - Z[z]);
    double result = y_0*(1-z_prop)+y_1*z_prop;
    return result;
}

/* Function to create a vector that is orthonormal to 
   the plane spanned by input vectors u and v */
void CrossProduct(const double *u, const double *v, double(&result)[3])
{
    result[0] = u[1]*v[2]-u[2]*v[1];
    result[1] = u[2]*v[0]-u[0]*v[2];
    result[2] = u[0]*v[1]-u[1]*v[0];
}

void WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData* NewImage(int width, int height)
{
    vtkImageData* image = vtkImageData::New();
    image->SetDimensions(PIXEL_WIDTH,PIXEL_HEIGHT,1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    return image;
}

/* Function to compute the dot product of two vectors */
double DotProduct(double* a, double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* The main procedure can be summarized as follows:
For each image (and its respective sample rate):
    For each pixel in the image:
        find the ray that goes from the camera through this pixel
        for each sample:
            intersect the volume with the ray for this pixel
            calculate color from intersection and previous samples
        assign resulting color to the pixel as output*/
int main(int argc, char* argv[])
{
    vtkDataSetReader* rdr = vtkDataSetReader::New();
    rdr->SetFileName("astro64.vtk");
    rdr->Update();
    time_t start_time = clock();
    int dims[3];
    vtkRectilinearGrid* rgrid = (vtkRectilinearGrid *)rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *)rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *)rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *)rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *)rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    struct Camera camera;
    camera = SetupCamera();
    struct ColorTransfer tf;
    tf = SetupColorTransfer();
    
    double magnitude, N[3], look[3], U[3], u[3], v[3], delta_x[3], delta_y[3];//we'll need several orthonormal vectors
    for (int i = 0; i < 3; i++)
        N[i] = camera.viewFocus[i] - camera.viewPosition[i];
    magnitude = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);//calculate length so we can normalize
    for (int i = 0; i < 3; i++)
        look[i] = N[i]/magnitude;//look is normalized n
    CrossProduct(look, camera.viewUp, U);//now U is orthonormal to our focus vector and orientation.
    magnitude = sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
    for (int i = 0; i < 3; i++)
        u[i] = U[i]/magnitude;//normalize U
    CrossProduct(look, u, v);//v is already normalized, as the cross product of two orthnormal vectors.
    for (int i = 0; i < 3; i++)
    {
        delta_x[i] = (((2*tan((camera.viewAngle/2) * (M_PI/180)))/PIXEL_WIDTH))*u[i];
        delta_y[i] = -(((2*tan((camera.viewAngle/2) * (M_PI/180)))/PIXEL_HEIGHT))*v[i];
    }
    
    vtkImageData *image;
    image = NewImage(PIXEL_WIDTH,PIXEL_HEIGHT);
    unsigned char *buffer;
    buffer = (unsigned char *)image->GetScalarPointer(0,0,0);
    for (int pic = 1; pic < NUM_PICS; pic++) {
        int num_samples = NUM_SAMPLES/(pic);
        printf("Working on picture number %d...\n",pic);
        printf("    Constructing ray for each pixel...\n");
        printf("        Sampling each ray %d times...\n", num_samples);
	    //for each pixel (i,j)
        for (int i = 0; i < PIXEL_WIDTH; i++)
        {
            for (int j = 0; j < PIXEL_HEIGHT; j++)
            {
                //construct the ray for this pixel
                double ray[3], position[3], sample_interval;
                for (int k = 0; k < 3; k++)
                    ray[k] = look[k] + ((((2*i)+1-PIXEL_WIDTH)/2)*delta_x[k]) + ((((2*j)+1-PIXEL_HEIGHT)/2)*delta_y[k]);
                magnitude = sqrt(ray[0]*ray[0]+ray[1]*ray[1]+ray[2]*ray[2]);//calculate length of ray so we can normalize
                for (int k = 0; k < 3; k++)
                    ray[k] = ray[k]/magnitude;//normalize the vector for pixel (i,j)
                for (int k = 0; k < 3; k++)
                    position[k] = camera.viewPosition[k] + ray[k]*camera.nearPlane;
                sample_interval = (camera.farPlane - camera.nearPlane)/(num_samples);//FIXME use length of hypotenuse if time is allowed (doesn't affect this data set)
                int first_sample_flag = 0;
                double running_color[3] = {0,0,0};
                double running_opacity = 0;
                
                //int check_flag = 0;
                for (int n = 0; n < (num_samples); n++) 
                {
                    for (int k = 0; k < 3; k++)
                        position[k] += ray[k]*sample_interval;
                    double s = EvaluateFieldAtLocation(dims,X,Y,Z,F,position);//FIXME brute force is slow
                    if (first_sample_flag == 0)//if this is the first colorable sample, initialize.
                    {
                        unsigned char *new_color = new unsigned char[3]();
                        first_sample_flag = 1;
                        tf.ApplyColorTransfer(s, new_color, running_opacity);
                        //apply alpha correction based on sample rate

                        running_opacity = 1.0 - pow((1.0 - running_opacity), (100.0/(num_samples)));
                        for (int idx = 0; idx < 3; idx++)
                            running_color[idx] = (double(new_color[idx]))*running_opacity;
                        free(new_color);
                    }
                    else
                    {
                        unsigned char *new_color = new unsigned char[3]();
                        double new_opacity = 0;
                        double alpha_to_go = 1 - running_opacity;
                        tf.ApplyColorTransfer(s, new_color, new_opacity);
                        //apply opacity correction to new alpha value
                        new_opacity = 1.0 - pow((1.0 - new_opacity), (100.0/(num_samples)));
                        for (int idx = 0; idx < 3; idx++)
                            running_color[idx] = running_color[idx] + alpha_to_go*new_opacity*(double(new_color[idx]));
                        running_opacity = running_opacity + alpha_to_go*new_opacity;
                        if (running_opacity >= 0.99)
                        {
                            break;
                        }
                        free(new_color);
                    }
                }//end for-each-sample loop

                int offset = 3*(j*PIXEL_WIDTH+i);
                unsigned char *place_in_buffer = buffer+offset;
                for (int k = 0; k < 3; k++)
                    place_in_buffer[k] = (unsigned char)running_color[k];
            }//end for each pixel j loop
        }//end for each pixel i loop
        
        cerr << "Execution time: " << (double)(clock() - start_time)/CLOCKS_PER_SEC << endl;
        char pic_num[40] = {0};
        char strbuf[10] = {0};
        strcat(pic_num, "star");
	if (pic < 10) {
		snprintf(strbuf, 10, "0%d", pic);
	} else { 
        	snprintf(strbuf, 10, "%d", pic);
        }
	strcat(pic_num,strbuf);
        WriteImage(image, pic_num);
	printf("Writing %s...\n", pic_num);
	start_time = clock();
    }
    return 0;
}
