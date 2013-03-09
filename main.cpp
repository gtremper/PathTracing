/* This is the main file for the raytracer */

#include <iostream>
#include <cstdio>
#include <stack>
#include <string>
#include <vector>
#include "omp.h"
#include <sstream>
#include <time.h>
#include <GLUT/glut.h>

#include "FreeImage.h"
#include "shaders.h"
#include "Shapes.h"
#include "Intersection.h"
#include "Light.h"
#include "Scene.h"

#define BPP 24
#define EPSILON 0.0000001
#define BUFFER_OFFSET(i) (reinterpret_cast<void*>(i))

using namespace std;

/* Paramaters */
double rays_cast;
Scene* scene;
int height;
int width;

/* Shaders */
GLuint vertexshader;
GLuint fragmentshader;
GLuint shaderprogram;
GLuint texture;
FIBITMAP* bitmap;
vec3 * pixels;
int update;

vec3 randHem(vec3& norm){
	double u1 = ((double)rand()/(double)RAND_MAX);
	double u2 = ((double)rand()/(double)RAND_MAX);
	
	double theta = acos(sqrt(u1));
	double phi = 2.0 * M_PI * u2;
	
	if (theta < EPSILON) {
		return norm;
	}
	if (1.0-norm[2] < EPSILON){
		return vec3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
	}
	
	const vec3 X = vec3(1,0,0);
	const vec3 Y = vec3(0,1,0);
	const vec3 Z = vec3(0,0,1);
	
	
	double dTheta = acos(glm::dot(norm, Z));
	vec3 xproj = glm::dot(X, norm) * X;
	vec3 yproj = glm::dot(Y,norm) * Y;
	vec3 xy_plane_vec = glm::normalize(xproj+yproj);
	double dPhi = acos(glm::dot(X, xy_plane_vec));
	if (xy_plane_vec[1] < 0.0){
		dPhi = -dPhi;
	}
	
	theta += dTheta;
	phi += dPhi;

	return vec3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));	
}

vec3 findColor(Scene* scene, Ray& ray, int depth) {

	Intersection hit = scene->KDTree->intersect(ray);
	
	if(!hit.primative || depth==1) {
		return vec3(0,0,0); //background color
	}
	
	vec3 normal = hit.primative->getNormal(hit.point);
	
	/* Tempararily return if hit light */
	if(hit.primative->emission[0]){
		return hit.primative->emission;
	}
	
	vec3 newDirection = randHem(normal);
	Ray newRay(hit.point+EPSILON*normal, newDirection);
	return hit.primative->diffuse * findColor(scene, newRay, depth-1);
	
	//double c1 = -glm::dot(normal, ray.direction);
	
	//vector<Light*>::iterator light = scene->lights.begin();
	//for(; light!=scene->lights.end(); ++light){
	//	color += (*light)->shade(hit, scene->KDTree, normal);
	//}
	
	//Ray reflectedRay = Ray(hit.point+EPSILON*normal, ray.direction+(2.0*normal*c1));
	
	//if(depth != 1) {
		//color += hit.primative->specular * findColor(scene, reflectedRay, depth-1);
		/*
		if(hit.primative->refractivity) {
			double n;
			if(c1 > 0.0) {
				n = 1.003f/hit.primative->indexofrefraction; // 1.009 is the refractive index of a vacuum
			} else {
				n = hit.primative->indexofrefraction/1.003f;
			}
			
			double c2 = sqrt(1 - n*n * (1 - c1*c1));
			Ray refractedRay = Ray(hit.point, glm::normalize((n*ray.direction) + (n*c1-c2)*normal));
			if(c1>0.0){
				refractedRay.origin -= EPSILON*normal;
			} else {
				refractedRay.direction = -refractedRay.direction;
				refractedRay.origin += EPSILON*normal;
			}
			color += hit.primative->refractivity * findColor(scene, refractedRay, depth-1);
		}
		*/
	//}
	//return color;
}

/* ouputs bitmap to global variable*/
void raytrace(double rayscast) {
	
	double subdivisions = 4;
	double subdivide = 1/subdivisions;
	
	double old_weight = rayscast/(rayscast+1.0);
	double new_weight = 1.0 - old_weight;
	
	#pragma omp parallel for
	for (int j=0; j<scene->height; j++){
		int tid = omp_get_thread_num();
		if(tid == 0) {
		   clog << "Progress: "<< (j*100*omp_get_num_threads())/scene->height <<"%"<<"\r";
		}
		RGBQUAD rgb;
		for (int i=0; i<scene->width; i++) {			
			vec3 color;
			for(double a=0; a<subdivisions; a+=1) {
				for(double b=0; b<subdivisions; b+=1) {
					double randomNum1 = ((double)rand()/(double)RAND_MAX) * subdivide;
					double randomNum2 = ((double)rand()/(double)RAND_MAX) * subdivide;
					Ray ray = scene->castEyeRay(i+(a*subdivide) + randomNum1,j+(b*subdivide)+randomNum2);
					//color += findColor(scene, ray, scene->maxdepth);
					color += findColor(scene, ray, 100);
				}
			}
			color *= (subdivide * subdivide);
			pixels[i + scene->width*j] *= old_weight;
			pixels[i + scene->width*j] += new_weight*color;
			color = pixels[i + scene->width*j];
			rgb.rgbRed = min(color[0],1.0)*255.0;
			rgb.rgbGreen = min(color[1],1.0)*255.0;
			rgb.rgbBlue = min(color[2],1.0)*255.0;
			FreeImage_SetPixelColor(bitmap,i,j,&rgb);
		}          
	}
}

void reshape(int w, int h){
	width = w;
	height = h;
	glViewport(0, 0, w, h);
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
	stringstream ss;
	switch(key){
		case 'l':
			update = 1;
			break;
		case 's':
			ss << scene->filename<<"_"<< int(rays_cast) << ".png";
			FreeImage_Save(FIF_PNG, bitmap, ss.str().c_str(), 0);
			cout << "Image saved!" << endl;
			break;
		case 'r':
			glutReshapeWindow(scene->width,scene->height);
			break;
		case 27:  // Escape to quit
			FreeImage_DeInitialise();
			delete scene;
			delete pixels;
			exit(0);
			break;
			
	}
	glutPostRedisplay();
}

void init(char* filename) {
	scene = new Scene(filename);
	rays_cast = 0.0;
	update = false;
	width = scene->width;
	height = scene->height;
	pixels = new vec3[width*height];
	memset(pixels, 0, sizeof(vec3)*width*height);
	FreeImage_Initialise();
	bitmap = FreeImage_Allocate(width, height, BPP);
	
	vertexshader = initshaders(GL_VERTEX_SHADER, "shaders/vert.glsl");
	fragmentshader = initshaders(GL_FRAGMENT_SHADER, "shaders/frag.glsl");
	shaderprogram = initprogram(vertexshader, fragmentshader);
	
	glGenTextures(1, &texture);
	glEnable(GL_TEXTURE_2D) ;
	glBindTexture(GL_TEXTURE_2D, texture);
	glActiveTexture(GL_TEXTURE0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	
	BYTE* bits = FreeImage_GetBits(bitmap);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, scene->width, scene->height,
		0, GL_BGR, GL_UNSIGNED_BYTE, (GLvoid*)bits);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1,1,-1,1,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glm::mat4 mv = glm::lookAt(glm::vec3(0,0,1),glm::vec3(0,0,0),glm::vec3(0,1,0));
	glLoadMatrixf(&mv[0][0]);
 
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT);
	
	if (update){
		time_t seconds = time(NULL);
		raytrace(rays_cast);
		BYTE* bits = FreeImage_GetBits(bitmap);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, scene->width, scene->height,
			0, GL_BGR, GL_UNSIGNED_BYTE, (GLvoid*)bits);
		rays_cast += 1.0;
		cout << "Number of Samples: " << rays_cast <<
		"\tTime: " << time(NULL)-seconds <<" seconds" << endl;
		update -= 1;
		glutPostRedisplay();
	}	
	
	glBegin(GL_QUADS);
	glTexCoord2d(0, 0); glVertex3d(-1, -1, 0);
	glTexCoord2d(0, 1); glVertex3d(-1, 1, 0);
	glTexCoord2d(1, 1); glVertex3d(1, 1, 0);
	glTexCoord2d(1, 0); glVertex3d(1, -1, 0);
	glEnd();
	
	glutSwapBuffers();
}

int main(int argc, char* argv[]){
	if(argc != 2) {
		cerr << "You need 1 scene file as the argument" << endl;
		exit(1);
	}
	srand(time(0));
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutCreateWindow("Path Tracer");
	init(argv[1]);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutReshapeWindow(width,height);
	glutMainLoop();
	return 0;
}