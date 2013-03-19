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
vec3 * direct_pixels;
int update;

const vec3 X = vec3(1,0,0);
const vec3 Y = vec3(0,1,0);
const vec3 Z = vec3(0,0,1);

/** rotate the Z vector in the direction of norm */
mat3 center_axis(const vec3& norm){
	vec3 perp = glm::cross(Z,norm);
	double angle = acos(glm::dot(Z,norm));
	mat4 rotation = glm::rotate(mat4(1.0), angle, perp);

	return mat3(rotation);
}

/* Sample a hemisphere for diffuse ray */
vec3 cos_weighted_hem(vec3& norm){
	double u1 = ((double)rand()/(double)RAND_MAX);
	double u2 = ((double)rand()/(double)RAND_MAX);

    vec3 y = vec3(norm);
    vec3 h = vec3(norm);
    double theta = acos(sqrt(1.0 - u1));
    double phi = 2.0 * M_PI * u2;
    double xs = sin(theta) * cos(phi);
    double ys = cos(theta);
    double zs = sin(theta) * sin(phi);
    if ((fabs(h[0]) <= fabs(h[1])) && (fabs(h[0]) <= fabs(h[2])))
      h[0] = 1.0;
    else if ((fabs(h[1]) <= fabs(h[0])) && (fabs(h[1]) <= fabs(h[2])))
      h[1] = 1.0;
    else
      h[2] = 1.0;
    vec3 x = glm::cross(h,y);
    vec3 z = glm::cross(x,y);

    vec3 direction = xs * x + ys * y + zs * z;
    return direction;

	//double theta = acos(sqrt(u1));
	//double phi = 2.0 * M_PI * u2;

	//if (theta < EPSILON) {
	//	return norm;
	//}

    //vec3 direction = vec3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));

    ///* Return direction rotated so its with respect to norm */
    //return center_axis(direction) * direction;
}

/* Sample a hemispehre for specular ray */
vec3 specular_weighted_hem(vec3& reflection, double n){
	double u1 = ((double)rand()/(double)RAND_MAX);
	double u2 = ((double)rand()/(double)RAND_MAX);

	double alpha = acos( pow( u1, 1.0 / (n + 1.0) ) );
	double phi = 2.0 * M_PI * u2;

	if (alpha < EPSILON) {
		return reflection;
	}
	vec3 direction = vec3(cos(phi)*sin(alpha), sin(phi)*sin(alpha), cos(alpha));

	/* return direction rotated so its with respecto to reflection */
	return center_axis(direction) * direction;
}

/* Shoot ray at scene */
vec3 findColor(Scene* scene, Ray& ray, int depth) {

	/* Intersect scene */
	Intersection hit = scene->KDTree->intersect(ray);

	// Should eventually replace this with russian roulette
	if(!hit.primative || depth==1) {
		return vec3(0,0,0); //background color
	}
	/* Tempararily return if hit light */
	// Maybe should weight this with specular and diffuse
	if( glm::length(hit.primative->emission) > EPSILON ){
		return hit.primative->emission;
	}

	vec3 normal = hit.primative->getNormal(hit.point);


	/*********************************************
	Importance sample global illumination
	*********************************************/
	double diffWeight = glm::length(hit.primative->diffuse);
	double specWeight = glm::length(hit.primative->specular);
	double threshold = diffWeight / (diffWeight + specWeight);
	double u1 = ((double)rand()/(double)RAND_MAX);

	vec3 color = vec3(0,0,0);
    //light_angle = (1.0 - (light_angle / (2 * M_PI)));
	light_angle = 1.0;
    Intersection newhit;
    vec3 newDirection;
    Ray newRay;
	/* Importance sample on macro level to choose diffuse or specular */
	if (u1 < threshold) {
		newDirection = cos_weighted_hem(normal);
		newRay = Ray(hit.point+EPSILON*normal, newDirection);
        newhit = scene->KDTree->intersect(newRay);
        if (newhit.primative) {
          if (glm::length(newhit.primative->emission) > EPSILON) {
            color = vec3(0,0,0); //ignore direct lighting contribution
          }
        }
		double prob = 1.0/threshold;
		color += prob * hit.primative->diffuse * findColor(scene, newRay, depth-1);
	} else {
		vec3 reflect = glm::reflect(ray.direction, normal);
		newDirection = specular_weighted_hem(reflect, hit.primative->shininess);
		newRay = Ray(hit.point+EPSILON*normal, newDirection);
        newhit = scene->KDTree->intersect(newRay);

		vec3 half = glm::normalize(hit.sourceDirection + newDirection);
		double phong =	pow( max(0.0,glm::dot(half,normal)) , hit.primative->shininess);

		/* Get probability for importance sampling */
		double cosalpha = glm::dot(newDirection,reflect);
		cosalpha = pow(cosalpha, hit.primative->shininess);
		double prob = cosalpha * (hit.primative->shininess + 1.0) / (2.0 * M_PI);

		double multiplier = phong / prob;
		multiplier *= 1.0/(1.0-threshold);
		color += multiplier * hit.primative->specular * max(0.0,glm::dot(normal, newDirection)) * findColor(scene, newRay, depth-1);
	}
	return color*M_PI;
}

/* Main raytracing function. Shoots ray for each pixel with anialiasing */
/* ouputs bitmap to global variable*/
void raytrace(double rayscast) {
	double subdivisions = scene->antialias;
	double subdivide = 1.0/subdivisions;

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
					color += findColor(scene, ray, scene->maxdepth);
				}
			}
			color *= (subdivide * subdivide);
			pixels[i + scene->width*j] *= old_weight;
			pixels[i + scene->width*j] += new_weight*color;
			color = pixels[i + scene->width*j] + direct_pixels[i + scene->width*j];
			rgb.rgbRed = min(color[0],1.0)*255.0;
			rgb.rgbGreen = min(color[1],1.0)*255.0;
			rgb.rgbBlue = min(color[2],1.0)*255.0;
			FreeImage_SetPixelColor(bitmap,i,j,&rgb);
		}
	}
}

/* Only the direct lighting */
vec3 direct_lighting(Scene* scene, Ray& ray){

	/* Intersect scene */
	Intersection hit = scene->KDTree->intersect(ray);

	if(!hit.primative) {
		return vec3(0,0,0); //background color
	}
	if( glm::length(hit.primative->emission) > EPSILON ){
		return hit.primative->emission;
	}

	vec3 color = vec3(0,0,0);
	for (unsigned int i = 0; i < scene->lights.size(); i++) {
		color += scene->lights[i]->shade(hit, scene->KDTree);
	}
	return color;
}


/* Calculate the direct lighting seperatly */
void direct_raytrace() {
	double subdivisions = scene->antialias;
	double subdivide = 1.0/subdivisions;

	#pragma omp parallel for
	for (int j=0; j<scene->height; j++){
		int tid = omp_get_thread_num();
		if(tid == 0) {
		   clog << "Direct Lighting: "<< (j*100*omp_get_num_threads())/scene->height <<"%"<<"\r";
		}
		RGBQUAD rgb;
		for (int i=0; i<scene->width; i++) {
			vec3 color;
			for(double a=0; a<subdivisions; a+=1) {
				for(double b=0; b<subdivisions; b+=1) {
					double randomNum1 = ((double)rand()/(double)RAND_MAX) * subdivide;
					double randomNum2 = ((double)rand()/(double)RAND_MAX) * subdivide;
					Ray ray = scene->castEyeRay(i+(a*subdivide) + randomNum1,j+(b*subdivide)+randomNum2);
					color += direct_lighting(scene, ray);
				}
			}
			color *= (subdivide * subdivide);
			direct_pixels[i + scene->width*j] = color;
			rgb.rgbRed = min(color[0],1.0)*255.0;
			rgb.rgbGreen = min(color[1],1.0)*255.0;
			rgb.rgbBlue = min(color[2],1.0)*255.0;
			FreeImage_SetPixelColor(bitmap,i,j,&rgb);
		}
	}
	clog << "	done\n";
}

/* Everything below here is openGL boilerplate */

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
		case 'f':
			update = 10000;
			break;
		case 'h':
			update = 0;
			break;
		case 'r':
			glutReshapeWindow(scene->width,scene->height);
			break;
		case 27:  // Escape to quit
			FreeImage_DeInitialise();
			delete scene;
			delete pixels;
			delete direct_pixels;
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
	direct_pixels = new vec3[width*height];
	memset(pixels, 0, sizeof(vec3)*width*height);
	memset(direct_pixels, 0, sizeof(vec3)*width*height);
	
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
		cout << "Iterations left: " << update << endl;
		time_t seconds = time(NULL);
		if (!rays_cast) {
			direct_raytrace();
		}
		//raytrace(rays_cast);
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
