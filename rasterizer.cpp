#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "hw2_types.h"
#include "hw2_math_ops.h"
#include "hw2_file_ops.h"
#include <iostream>
#include <math.h>
#include <algorithm>
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

typedef struct
{
	int x, y, colorId;
} Vec2;

Color addColor(const Color& a, const Color& b) {
	Color result;
	result.r = a.r + b.r;
	result.g = a.g + b.g;
	result.b = a.b + b.b;
	return result;
}
Color subtractColor(const Color& a, const Color& b) {
	Color result;
	result.r = a.r - b.r;
	result.g = a.g - b.g;
	result.b = a.b - b.b;
	return result;
}

Color divideWScalar(const Color& c, double s) {
	Color result;
	result.r = c.r / s;
	result.g = c.g / s;
	result.b = c.b / s;
	return result;
}
Color multiplyWScalar(const Color& c, double s) {
	Color result;
	result.r = c.r * s;
	result.g = c.g * s;
	result.b = c.b * s;
	return result;
}

void rotateX(double theta, double m[4][4]) {
	makeIdentityMatrix(m);
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
	makeIdentityMatrix(m);

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
	makeIdentityMatrix(m);

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
	double tmpResult[4][4];
	multiplyMatrixWithMatrix(tmpResult, m1, m2);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			m1[i][j] = tmpResult[i][j];
		}
	}
}

void rotateArbitrary(double theta, double a, double b, double c, double r[4][4]) { // to be tested
	if ((a == 1 && b == 0) && c == 0) {
		rotateX(theta, r);
	} else if ((a == 0 && b == 1) && c == 0) {
		rotateY(theta, r);
	} else if ((a == 0 && b == 0) && c == 1) {
		rotateZ(theta, r);
	} else {
		double rotateXMatrix[4][4];
		double rotateYMatrix[4][4];
		double rotateZMatrix[4][4];
		double rotateYMatrixM[4][4];
		double rotateXMatrixM[4][4];
		double bcs = sqrt(b*b + c * c);
		double alpha = acos(c / (bcs))* 180.0 / PI;
		double beta = acos(bcs)* 180.0 / PI;
		// Rotations
		rotateX(-alpha, rotateXMatrix);
		rotateY(beta, rotateYMatrix);
		rotateZ(theta, rotateZMatrix);
		rotateY(-beta, rotateYMatrixM);
		rotateX(alpha, rotateXMatrixM);
		// Create composite matrix
		assignedMultiplication(r, rotateXMatrix);
		assignedMultiplication(r, rotateYMatrix);
		assignedMultiplication(r, rotateZMatrix);
		assignedMultiplication(r, rotateYMatrixM);
		assignedMultiplication(r, rotateXMatrixM);
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

void createModelTransformationMatrix(Model& m, double transMatrix[4][4]) {
	makeIdentityMatrix(transMatrix);
	for (int i = m.numberOfTransformations - 1; i >= 0; i--) {
		double tmpMatrix[4][4];
		makeIdentityMatrix(tmpMatrix);
		if (m.transformationTypes[i] == 'r') {
			auto& rotation = rotations[m.transformationIDs[i]];
			rotateArbitrary(rotation.angle, rotation.ux, rotation.uy, rotation.uz, tmpMatrix);

		} else if (m.transformationTypes[i] == 't') {
			auto translation = translations[m.transformationIDs[i]];
			translate(translation, tmpMatrix);

		} else if (m.transformationTypes[i] == 's') {
			auto scaling = scalings[m.transformationIDs[i]];
			scale(scaling, tmpMatrix);
		}
		assignedMultiplication(transMatrix, tmpMatrix);
	}
}

void createCameraTransformationMatrix(Camera& cam, double cameraMatrix[4][4]) {
	cameraMatrix[0][0] = cam.u.x;
	cameraMatrix[0][1] = cam.u.y;
	cameraMatrix[0][2] = cam.u.z;
	cameraMatrix[0][3] = -(cam.u.x * cam.pos.x + cam.u.y*cam.pos.y + cam.u.z* cam.pos.z);
	cameraMatrix[1][0] = cam.v.x;
	cameraMatrix[1][1] = cam.v.y;
	cameraMatrix[1][2] = cam.v.z;
	cameraMatrix[1][3] = -(cam.v.x * cam.pos.x + cam.v.y*cam.pos.y + cam.v.z* cam.pos.z);
	cameraMatrix[2][0] = cam.w.x;
	cameraMatrix[2][1] = cam.w.y;
	cameraMatrix[2][2] = cam.w.z;
	cameraMatrix[2][3] = -(cam.w.x * cam.pos.x + cam.w.y*cam.pos.y + cam.w.z* cam.pos.z);
	cameraMatrix[3][0] = 0;
	cameraMatrix[3][1] = 0;
	cameraMatrix[3][2] = 0;
	cameraMatrix[3][3] = 1;
}

void perspectiveDivide(double d[4]) {
	d[0] /= d[3];
	d[1] /= d[3];
	d[2] /= d[3];
	d[3] /= d[3];
}

// Takes two vertices and draws the line between them NOTE:VERTEX2 IS NE OF VERTEX1
void midPointE_NE(const Vec2& vertex1, const Vec2& vertex2) {
	int y, d, dInc1, dInc2;
	Color c, cInc;
	y = vertex1.y;
	d = 2 * (vertex1.y - vertex2.y) + (vertex2.x - vertex1.x);
	c = colors[vertex1.colorId];
	dInc1 = 2 * (vertex1.y - vertex2.y) + 2 * (vertex2.x - vertex1.x);
	dInc2 = 2 * (vertex1.y - vertex2.y);
	cInc = divideWScalar(subtractColor(colors[vertex2.colorId],c), vertex2.x - vertex1.x);
	for (int x = vertex1.x; x < vertex2.x; x++) {
		image[x][y] = c;
		if (d < 0) { // chose NE
			y++;
			d += dInc1;
		} else { // chose E
			d += dInc2;
		}
		c = addColor(c, cInc);
	}
}

// Takes two vertices and draws the line between them NOTE:VERTEX2 IS NE OF VERTEX1 (TO BE CHECKED)
void midPointN_NE(const Vec2& vertex1, const Vec2& vertex2) {

	int x, d, dInc1, dInc2;
	Color c, cInc;
	x = vertex1.x;
	d = (vertex2.y - vertex1.y) + 2 * (vertex1.x - vertex2.x);
	c = colors[vertex1.colorId];
	dInc1 = 2 * (vertex2.y - vertex1.y) + 2 * (vertex1.x - vertex2.x);
	dInc2 = 2 * (vertex1.x - vertex2.x);
	cInc = divideWScalar(subtractColor(colors[vertex2.colorId], c), vertex2.y - vertex1.y);
	for (int y = vertex1.y; y < vertex2.y; y++) {
		image[x][y] = c;
		if (d < 0) { // chose NE

			x++;

			d += dInc1;
		} else { // chose N
			d += dInc2;
		}
		c = addColor(c, cInc);
	}
}

// Takes two vertices and draws the line between them NOTE:VERTEX2 IS NW OF VERTEX1 (TO BE CHECKED)
void midPointW_NW(const Vec2& vertex1, const Vec2& vertex2) {
	//std::cout << "Test v1:" << vertex1.x << "  " << vertex1.y << " v2:" << vertex2.x << " " << vertex2.y << std::endl;

	int y, d, dInc1, dInc2;
	Color c, cInc;
	y = vertex1.y;
	d = 2 * (vertex2.y - vertex1.y) - (vertex1.x - vertex2.x);
	c = colors[vertex1.colorId];
	dInc1 = 2 * (vertex2.y - vertex1.y) - 2 * (vertex1.x - vertex2.x);
	dInc2 = 2 * (vertex2.y - vertex1.y);
	cInc = divideWScalar(subtractColor(colors[vertex2.colorId], c), vertex2.x - vertex1.x);
	for (int x = vertex1.x; x < vertex2.x; x++) {
		image[x][y] = c;
		if (d < 0) { // chose NW
			y--;
			d += dInc1;
		} else { // chose W
			d += dInc2;
		}
		c = addColor(c, cInc);
	}
}

// Takes two vertices and draws the line between them NOTE:VERTEX2 IS NW OF VERTEX1 (TO BE CHECKED)
void midPointN_NW(const Vec2& vertex1, const Vec2& vertex2) {

	int x, d, dInc1, dInc2;
	Color c, cInc;
	x = vertex1.x;
	d = 2 * (vertex2.x - vertex1.x) - (vertex1.y - vertex2.y);
	c = colors[vertex1.colorId];
	dInc1 = 2 * (vertex2.x - vertex1.x) - 2 * (vertex1.y - vertex2.y);
	dInc2 = 2 * (vertex2.x - vertex1.x);
	cInc = divideWScalar(subtractColor(colors[vertex2.colorId], c), vertex2.y - vertex1.y);
	for (int y = vertex1.y; y < vertex2.y; y++) {
		image[x][y] = c;
		if (d < 0) { // chose NW
			x--;
			d += dInc1;
		} else { // chose N
			d += dInc2;
		}
		c = addColor(c, cInc);
	}
}

// Takes two vertices and draws the line between them NOTE:VERTEX2 IS NORTH OF VERTEX1
void vertical(const Vec2& vertex1, const Vec2& vertex2) {
	Color c, cInc;
	c = colors[vertex1.colorId];
	cInc = divideWScalar(subtractColor(colors[vertex2.colorId], c), vertex2.y - vertex1.y);

	for (int y = vertex1.y; y < vertex2.y; y++) {
		image[vertex1.x][y] = c;
		c = addColor(c, cInc);
	}
}

// Takes two vertices and draws the line between them NOTE:VERTEX2 IS EAST OF VERTEX1
void horizontal(const Vec2& vertex1, const Vec2& vertex2) {
	Color c, cInc;
	c = colors[vertex1.colorId];
	cInc = divideWScalar(subtractColor(colors[vertex2.colorId], c), vertex2.x - vertex1.x);

	for (int x = vertex1.x; x < vertex2.x; x++) {
		image[x][vertex1.y] = c;
		c = addColor(c, cInc);
	}
}

// Takes two arbitrary vertices and draws the line between them
void lineRasterizer(const Vec2& vertex1, const Vec2& vertex2) {

	int slopeNom = vertex2.y - vertex1.y, slopeDenom = vertex2.x - vertex1.x;
	if (slopeDenom == 0) {
		if (slopeNom > 0) vertical(vertex1, vertex2);
		else vertical(vertex2, vertex1);
		return;
	} else if (slopeNom == 0) {
		if (slopeDenom > 0) horizontal(vertex1, vertex2);
		else horizontal(vertex2, vertex1);
		return;
	}
	auto slope = ((float)slopeNom) / (slopeDenom);
	if (slope > 1) {
		if (slopeNom > 0) midPointN_NE(vertex1, vertex2);
		else midPointN_NE(vertex2, vertex1);
	} else if (0 < slope && slope < 1) {
		if (slopeDenom > 0) midPointE_NE(vertex1, vertex2);
		else midPointE_NE(vertex2, vertex1);
	} else if (0 > slope && slope > -1) {
		if (slopeDenom > 0) midPointW_NW(vertex1, vertex2);
		else midPointW_NW(vertex2, vertex1);
	} else if (slope < -1) {
		if (slopeNom > 0) midPointN_NW(vertex1, vertex2);
		else midPointN_NW(vertex2, vertex1);
	}

}
inline double f(double x, double y, const Vec2& v_0, const  Vec2& v_1) noexcept {
	return x * (v_0.y - v_1.y) + y * (v_1.x - v_0.x) + v_0.x*v_1.y - v_0.y* v_1.x;
}
inline Vec2 convertToVec2(double vp[4], Vec3& v) noexcept {
	return Vec2{ (int)vp[0], (int)vp[1], v.colorId };
}
void triangleRasterize(const Vec2& v0, const Vec2& v1, const Vec2& v2) {
	auto y_min = v0.y < v1.y ? (v0.y < v2.y ? v0.y : v2.y) : (v1.y < v2.y ? v1.y : v2.y);
	auto x_min = v0.x < v1.x ? (v0.x < v2.x ? v0.x : v2.x) : (v1.x < v2.x ? v1.x : v2.x);
	auto y_max = v0.y > v1.y ? (v0.y > v2.y ? v0.y : v2.y) : (v1.y > v2.y ? v1.y : v2.y);
	auto x_max = v0.x > v1.x ? (v0.x > v2.x ? v0.x : v2.x) : (v1.x > v2.x ? v1.x : v2.x);
	for (int y = y_min; y < y_max; y++) {
		for (int x = x_min; x < x_max; x++) {
		
			auto alpha = f(x, y, v1, v2) / f(v0.x, v0.y, v1, v2);
			auto beta = f(x, y, v2, v0) / f(v1.x, v1.y, v2, v0);
			auto gamma = f(x, y, v0, v1) / f(v2.x, v2.y, v0, v1);
			if ((alpha >= 0 && beta >= 0) && gamma >= 0) {
				auto ca = multiplyWScalar(colors[v0.colorId], alpha);
				auto cb = multiplyWScalar(colors[v1.colorId], beta);
				auto cg = multiplyWScalar(colors[v2.colorId], gamma);
				auto colorSum = addColor(addColor(ca, cb), cg);

				image[x][y] = colorSum;
			}
		}
	}
}
inline Vec3 getTriangleNormal(Vec3 u, Vec3 v) noexcept {
	return crossProductVec3(u,v);
}
bool solidCulling(Vec3& v, Vec3&n) {
	return (v.x*n.x + v.y*n.y + v.z* n.z) > 0;
}
void forwardRenderingPipeline(Camera cam) {

	for (int i = 0; i <= numberOfModels; i++) {
		auto & model = models[i];
		double resultingMatrix[4][4];
		makeIdentityMatrix(resultingMatrix);
		double cameraCoord[4][4];
		makeIdentityMatrix(cameraCoord);
		double cameraMatrix[4][4];
		double transMatrix[4][4];
		double perspectiveMatrix[4][4] = { {2 * abs(cam.n) / (cam.r - cam.l),0 , (cam.r + cam.l) / (cam.r - cam.l),0},
										   {0, 2 * abs(cam.n) / (cam.t - cam.b), (cam.t + cam.b) / (cam.t - cam.b), 0},
										   {0,0, (abs(cam.n) + abs(cam.f)) / (abs(cam.n) - abs(cam.f)), 2 * abs(cam.f)*abs(cam.n) / (abs(cam.n) - abs(cam.f))},
										   {0,0,-1,0} };
		double viewPortMatrix[4][4] = { {cam.sizeX / 2.0,0,0.0, (cam.sizeX - 1) / 2.0},
										 {0.0, cam.sizeY / 2.0,0 , (cam.sizeY - 1) / 2.0},
										 {0.0,0.0,1.0,0.0},
										 {0.0,0.0,0.0,1.0} };
		createModelTransformationMatrix(model, transMatrix);
		createCameraTransformationMatrix(cam, cameraMatrix);
		// Compose the matrices
		assignedMultiplication(cameraCoord, perspectiveMatrix);
		assignedMultiplication(cameraCoord, cameraMatrix);
		assignedMultiplication(cameraCoord, transMatrix);
		assignedMultiplication(resultingMatrix, viewPortMatrix);
		assignedMultiplication(resultingMatrix, perspectiveMatrix);
		assignedMultiplication(resultingMatrix, cameraMatrix);
		assignedMultiplication(resultingMatrix, transMatrix);
		std::vector<Vec2> vec;
		std::vector<Vec3> cullVec;
		for (int i = 0; i < model.numberOfTriangles; i++) {
			auto &modelVertices = model.triangles[i];
			if (backfaceCullingSetting) {
				for (int j = 0; j < 3; j++) {
					auto &vertex = vertices[modelVertices.vertexIds[j]];
					double v[4] = { vertex.x, vertex.y, vertex.z, 1 };
					double r[4];
					multiplyMatrixWithVec4d(r, cameraCoord, v);
					cullVec.emplace_back(Vec3{ r[0],r[1],r[2] });

				}
				auto normal = normalizeVec3(getTriangleNormal(subtractVec3(cullVec[2], cullVec[0]), subtractVec3(cullVec[1], cullVec[0])));
				auto center = Vec3{ (cullVec[0].x + cullVec[1].x + cullVec[2].x) / 3.0, (cullVec[0].y + cullVec[1].y + cullVec[2].y) / 3.0 ,(cullVec[0].z + cullVec[1].z + cullVec[2].z) / 3.0 };
				if (solidCulling(center, normal)) {
					cullVec.clear();
					continue;
				}
			}

			for (int j = 0; j < 3; j++) {
				auto &vertex = vertices[modelVertices.vertexIds[j]];
				double v[4] = { vertex.x, vertex.y, vertex.z, 1 };
				double r[4];
				multiplyMatrixWithVec4d(r, resultingMatrix, v);
				perspectiveDivide(r);
				vec.emplace_back(convertToVec2(r, vertex));

			}

			if (model.type) {
				triangleRasterize(vec[0], vec[1], vec[2]);
			} 
			else {
				lineRasterizer(vec[0], vec[1]);
				lineRasterizer(vec[1], vec[2]);
				lineRasterizer(vec[2], vec[0]);
			}

			vec.clear();
			cullVec.clear();

		}

		//drawline with perspective divide
		// Continue with 'resultingMatrix'
	}
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
