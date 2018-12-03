#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "hw2_types.h"
#include "hw2_math_ops.h"
#include "hw2_file_ops.h"
#include <iostream>
#include <math.h>
#include <vector>
#define PI 3.14159265
Camera cameras[100];
int numberOfCameras = 0;

Model models[1000];
int numberOfModels = 0;

Color colors[100000];
int numberOfColors = 0;

Translation translations[1000];
int numberOfTranslations = 0;

Rotation rotations[1000];
int numberOfRotations = 0;

Scaling scalings[1000];
int numberOfScalings = 0;

Vec3 vertices[100000];
int numberOfVertices = 0;

Color backgroundColor;

// backface culling setting, default disabled
int backfaceCullingSetting = 0;

Color **image;



/*
	Initializes image with background color
*/
void initializeImage(Camera cam) {
	int i, j;

	for (i = 0; i < cam.sizeX; i++)
		for (j = 0; j < cam.sizeY; j++) {
			image[i][j].r = backgroundColor.r;
			image[i][j].g = backgroundColor.g;
			image[i][j].b = backgroundColor.b;

		}
}

/*
	Transformations, culling, rasterization are done here.
	You can define helper functions inside this file (rasterizer.cpp) only.
	Using types in "hw2_types.h" and functions in "hw2_math_ops.cpp" will speed you up while working.
*/
void rotateX(double theta,double m[4][4]) {
	double cosVal = cos(theta * PI / 180.0);
	double sinVal = sin(theta * PI / 180.0);
	m[0][0] = 1;
	m[1][1] = cosVal;
	m[1][2] = -sinVal;
	m[2][1] = sinVal;
	m[2][2] = cosVal;
	m[3][3] = 1;
}
void rotateY(double theta, double m[4][4]) {
	double cosVal = cos(theta * PI / 180.0);
	double sinVal = sin(theta * PI / 180.0);
	m[0][0] = cosVal;
	m[0][2] = sinVal;
	m[1][1] = 1;
	m[2][0] = -sinVal;
	m[2][2] = cosVal;
	m[3][3] = 1;
}
void rotateZ(double theta, double m[4][4]) {
	double cosVal = cos(theta * PI / 180.0);
	double sinVal = sin(theta * PI / 180.0);
	m[0][0] = cosVal;
	m[0][1] = -sinVal;
	m[1][0] = sinVal;
	m[1][1] = cosVal;
	m[2][2] = 1;
	m[3][3] = 1;
}
void assignedMultiplication(double m1[4][4], double m2[4][4]) {
	double tmpResult[4][4] = { {0} };
	multiplyMatrixWithMatrix(tmpResult, m1, m2);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			m1[i][j] = tmpResult[i][j];
		}
	}
}
void rotateArbitrary(double theta, double a, double b ,double c, double r[4][4]) {
	if ((a == 1 && b == 0) && c == 0) {
		rotateX(theta, r);
	} else if ((a == 0 && b == 1) && c == 0) {
		rotateY(theta, r);
	} else if ((a == 0 && b == 0) && c == 1) {
		rotateZ(theta, r);
	} else {
		double rotateXMatrix[4][4] = { {0} };
		double rotateYMatrix[4][4] = { {0} };
		double rotateZMatrix[4][4] = { {0} };
		double res1[4][4] = { {0} };
		double res2[4][4] = { {0} };
		double res3[4][4] = { {0} };
		double bcs = sqrt(b*b + c * c);
		double alpha = acos(c / (bcs))* 180.0 / PI;
		double beta = acos(bcs)* 180.0 / PI;
		rotateX(alpha, rotateXMatrix);
		rotateY(-beta, rotateYMatrix);
		rotateZ(theta, rotateZMatrix);
		multiplyMatrixWithMatrix(res1, rotateXMatrix, rotateYMatrix);
		multiplyMatrixWithMatrix(res2, res1, rotateZMatrix);
		rotateY(beta, rotateYMatrix);
		rotateX(-alpha, rotateXMatrix);
		multiplyMatrixWithMatrix(res3, res2, rotateYMatrix);
		multiplyMatrixWithMatrix(r, res3, rotateXMatrix);
	}
	
}
void translate(Translation& t, double transMatrix[4][4]) {
	transMatrix[0][3] = t.tx;
	transMatrix[1][3] = t.ty;
	transMatrix[2][3] = t.tz;
}
void scale(Scaling&s, double scaleMatrix[4][4]) {
	scaleMatrix[0][0] = s.sx;
	scaleMatrix[1][1] = s.sy;
	scaleMatrix[2][2] = s.sz;
}

void createModelTransformationMatrix(Model& m , double transMatrix[4][4]) {
	makeIdentityMatrix(transMatrix);
	for (int i = 0; i < m.numberOfTransformations; i++) {
		if (m.transformationTypes[i] == 'r') {
			double tmpMatrix[4][4];
			makeIdentityMatrix(tmpMatrix);
			auto rotation = rotations[m.transformationIDs[i]];
			rotateArbitrary(rotation.angle, rotation.ux, rotation.uy, rotation.uz, tmpMatrix);
			assignedMultiplication(transMatrix, tmpMatrix);
		} else if (m.transformationTypes[i] == 't') {
			double tmpMatrix[4][4];
			makeIdentityMatrix(tmpMatrix);
			auto translation = translations[m.transformationIDs[i]];
			translate(translation, tmpMatrix);
			assignedMultiplication(transMatrix, tmpMatrix);
		} else if (m.transformationTypes[i] == 's') {
			double tmpMatrix[4][4];
			makeIdentityMatrix(tmpMatrix);
			auto scaling = scalings[m.transformationIDs[i]];
			scale(scaling, tmpMatrix);
			assignedMultiplication(transMatrix, tmpMatrix);
		}
	}
	
}

void createCameraTransformationMatrix(Camera& cam, double cameraMatrix[4][4]) {
	double translationMatrix[4][4];
	makeIdentityMatrix(translationMatrix);
	makeIdentityMatrix(cameraMatrix);
	// translate to the origin
	Translation t{-cam.pos.x,-cam.pos.y,-cam.pos.z};
	translate(t, translationMatrix);
	assignedMultiplication(cameraMatrix, translationMatrix);
	double rotationMatrix[4][4] = { {cam.u.x,cam.u.y, cam.u.z, 0},{cam.v.x,cam.v.y,cam.v.z,0},
                                    {cam.w.x,cam.w.y,cam.w.z,0},{0,0,0,1} };
	assignedMultiplication(cameraMatrix, rotationMatrix);
}
void perspectiveDivide(double d[4]) {
	d[0] /= d[3];
	d[1] /= d[3];
	d[2] /= d[3];
	d[3] /= d[3];
}
void forwardRenderingPipeline(Camera cam) {
	for (int i = 0; i <= numberOfModels; i++) {
		double resultingMatrix[4][4];
		makeIdentityMatrix(resultingMatrix);
		double transMatrix[4][4];
		double cameraMatrix[4][4];
		double perspectiveMatrix[4][4] = { {2 * cam.n / (cam.r - cam.l),0 , (cam.r + cam.l) / (cam.r - cam.l),0},
										   {0, 2 * cam.n / (cam.t - cam.b), (cam.t + cam.b) / (cam.t - cam.b), 0},
										   {0,0,-(cam.f + cam.n) / (cam.f - cam.n), -2 * cam.f*cam.n / (cam.f - cam.n)},
										   {0,0,-1,0}};
		createModelTransformationMatrix(models[i],transMatrix);
		createCameraTransformationMatrix(cam, cameraMatrix);
		// Multiply with 4d vector
		assignedMultiplication(resultingMatrix, transMatrix);
		assignedMultiplication(resultingMatrix, cameraMatrix);
		assignedMultiplication(resultingMatrix, perspectiveMatrix);
		// Perspective divide multiplication
		// Viewport transformation
		int b = 0;
	}
	
	// TODO: IMPLEMENT HERE
}


int main(int argc, char **argv) {
	int i, j;

	if (argc < 2) {
		std::cout << "Usage: ./rasterizer <scene file> <camera file>" << std::endl;
		return 1;
	}

	// read camera and scene files
	readSceneFile(argv[1]);
	readCameraFile(argv[2]);

	image = 0;

	for (i = 0; i < numberOfCameras; i++) {

		// allocate memory for image
		if (image) {
			for (j = 0; j < cameras[i].sizeX; j++) {
				delete image[j];
			}

			delete[] image;
		}

		image = new Color*[cameras[i].sizeX];

		if (image == NULL) {
			std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
			exit(1);
		}

		for (j = 0; j < cameras[i].sizeX; j++) {
			image[j] = new Color[cameras[i].sizeY];
			if (image[j] == NULL) {
				std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
				exit(1);
			}
		}


		// initialize image with basic values
		initializeImage(cameras[i]);

		// do forward rendering pipeline operations
		forwardRenderingPipeline(cameras[i]);

		// generate PPM file
		writeImageToPPMFile(cameras[i]);

		// Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
		// Notice that os_type is not given as 1 (Ubuntu) or 2 (Windows), below call doesn't do conversion.
		// Change os_type to 1 or 2, after being sure that you have ImageMagick installed.
		convertPPMToPNG(cameras[i].outputFileName, 99);
	}
	system("pause");
	return 0;

}
