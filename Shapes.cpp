/* Class definitions for Triangles and Spheres */
#include <iostream>
#include <vector>

#include "Shapes.h"
#include "Intersection.h"

#define EPSILON 0.000001f

using namespace std;

/***  RAY  ***/
Ray::Ray(const vec3& o,const vec3& d){
	origin = o;
	direction = glm::normalize(d);
}

vec3 Ray::getPoint(float t) {
	return origin + t*direction;
}

/***  INTERSECTION  ***/
Intersection::Intersection(vector<Shape*>& objects, Ray& ray) {
	float min_t = DBL_MAX;
	primative = NULL;
	for(vector<Shape*>::iterator prim=objects.begin(); prim!=objects.end(); prim++) {
		float t = (*prim)->intersect(ray);
		if (t>0 && t<min_t){
			primative = *prim;
			min_t = t;
		}
	}
	point = ray.getPoint(min_t);
	sourceDirection = -ray.direction;
}

/***  TRIANGLE  ***/
Triangle::Triangle(vec3 point0, vec3 point1, vec3 point2) {
	p0 = point0;
	p1 = point1;
	p2 = point2;
	n0 = glm::normalize(glm::cross(p1-p0,p2-p0));
	float minx = min(min(p0[0],p1[0]),p2[0]);
	float maxx = max(max(p0[0],p1[0]),p2[0]);
	float miny = min(min(p0[1],p1[1]),p2[1]);
	float maxy = max(max(p0[1],p1[1]),p2[1]);
	float minz = min(min(p0[2],p1[2]),p2[2]);
	float maxz = max(max(p0[2],p1[2]),p2[2]);
	vec3 minvec = vec3(minx,miny,minz);
	vec3 maxvec = vec3(maxx,maxy,maxz);
	aabb = AABB(minvec, maxvec);
}


float Triangle::intersect(Ray& ray){
	vec3 col1 = p1-p0;
	vec3 col2 = p2-p0;
	mat3 M = mat3(col1, col2, -ray.direction);
	float det = glm::determinant(M);
	if (det<EPSILON && det>-EPSILON) return -1.0f;
	
	M[0] = ray.origin-p0;
	float alpha = glm::determinant(M)/det;
	if (alpha<0.0f || alpha>1.0f) return -1.0f;
	M[0] = col1;
	M[1] = ray.origin-p0;

	float beta = glm::determinant(M)/det;
	if (beta<0.0f || beta+alpha>1.0f) return -1.0f;
	M[1] = col2;
	M[2] = ray.origin-p0;
	return glm::determinant(M)/det;
}

vec3 Triangle::getNormal(vec3& hit){
	return n0;
}
/*
vec3 Triangle::getTexture(vec3& hit){
	mat2 M = mat2(p1[0]-p0[0], p1[1]-p0[1], p2[0]-p0[0], p2[1]-p0[1]);
	float det = glm::determinant(M);
	M[0][0] = hit[0]-p0[0];
	M[0][1] = hit[1]-p0[1];
	float alpha = glm::determinant(M)/det;//p0 to p1
	M[0][0] = p1[0]-p0[0];
	M[0][1] = p1[1]-p0[1];
	M[1][0] = hit[0]-p0[0];
	M[1][1] = hit[1]-p0[1];
	float beta = glm::determinant(M)/det;//p0 to p2
}

*/
/***  NORMTRIANGLE  ***/
NormTriangle::NormTriangle(vec3 point0, vec3 point1, vec3 point2,
							vec3 norm0, vec3 norm1, vec3 norm2) {
	p0 = point0;
	p1 = point1;
	p2 = point2;
	n0 = glm::normalize(norm0);
	n1 = glm::normalize(norm1);
	n2 = glm::normalize(norm2);
}

vec3 NormTriangle::getNormal(vec3& hit){
	mat2 M = mat2(p1[0]-p0[0], p1[1]-p0[1], p2[0]-p0[0], p2[1]-p0[1]);
	float det = glm::determinant(M);
	M[0][0] = hit[0]-p0[0];
	M[0][1] = hit[1]-p0[1];
	float beta = glm::determinant(M)/det;
	M[0][0] = p1[0]-p0[0];
	M[0][1] = p1[1]-p0[1];
	M[1][0] = hit[0]-p0[0];
	M[1][1] = hit[1]-p0[1];
	float gamma = glm::determinant(M)/det;
	return glm::normalize((1-beta-gamma)*n0 + beta*n1 + gamma*n2);
}

/***  SPHERE  ***/
Sphere::Sphere(mat4 trans){
	mv = trans;
	inv = glm::inverse(mv);
	
	/*Bounding Box*/
	mat4 S = mat4(1.0f);
	S[3][3] = -1.0f;
	mat4 R = mv*S*glm::transpose(mv);
	float minx = (R[0][3]+sqrt(R[0][3]*R[0][3] -R[3][3]*R[0][0]))/R[3][3];
	float maxx = (R[0][3]-sqrt(R[0][3]*R[0][3] -R[3][3]*R[0][0]))/R[3][3];
	float miny = (R[1][3]+sqrt(R[1][3]*R[1][3] -R[3][3]*R[1][1]))/R[3][3];
	float maxy = (R[1][3]-sqrt(R[1][3]*R[1][3] -R[3][3]*R[1][1]))/R[3][3];
	float minz = (R[2][3]+sqrt(R[2][3]*R[2][3] -R[3][3]*R[2][2]))/R[3][3];
	float maxz = (R[2][3]-sqrt(R[2][3]*R[2][3] -R[3][3]*R[2][2]))/R[3][3];
	
	vec3 minvec = vec3(minx,miny,minz);
	vec3 maxvec = vec3(maxx,maxy,maxz);
	aabb = AABB(minvec, maxvec);
}



float Sphere::intersect(Ray& ray) {
	vec3 direction = glm::normalize(vec3(inv * vec4(ray.direction,0.0f)));
	vec3 origin =vec3(inv * vec4(ray.origin,1.0f));
	
	float b = 2.0f * glm::dot(direction, origin);
	float c = glm::dot(origin,origin) - 1.0f;
	float det = b*b - 4.0f*c;
	if (det<EPSILON) return -1.0f;
	det = sqrt(det);
	
	float t1 = (-b+det)*0.5f;
	float t2 = (-b-det)*0.5f;
	
	if (t1<0.0 && t2<0.0) return -1.0f;
	
	if (t2>0.0f) {
		vec4 hit = mv * vec4(origin+t2*direction,1.0f);
		return glm::distance(ray.origin,vec3(hit));
	} else { //t1 is closer
		vec4 hit = mv * vec4(origin+t1*direction,1.0f);
		return glm::distance(ray.origin,vec3(hit));
	}
}

vec3 Sphere::getNormal(vec3& hit){
	return glm::normalize(vec3(glm::transpose(inv)*inv*vec4(hit,1.0f)));
}

/***  AABB  ***/
AABB::AABB(vec3& minarg, vec3& maxarg) {
	this->aabbmin = minarg;
	this->aabbmax = maxarg;
	this->center = (minarg+maxarg)/2.0f;
	
}

inline
bool intersect1D(float start, float dir, float axisMin, float axisMax, float& near, float& far){
	// Parallel
	if(dir<EPSILON && dir>-EPSILON){
		return (start>axisMin) && (start<axisMax);
	}
	
	//intersection parameters
	float t0 = (axisMin-start)/dir;
	float t1 = (axisMax-start)/dir;
	
	if(t0>t1){
		float temp = t1;
		t1 = t0;
		t0 = temp;
	}
	
	near = max(t0,near);
	far = min(t1,far);
	
	if(near>far || far<0.0f) return false;
	
	return true;
}

bool
operator<(const vec3 &vecA, const vec3 &vecB){ 
	if (vecA[0] > vecB[0]) return false;
	if (vecA[1] > vecB[1]) return false;
	if (vecA[2] > vecB[2]) return false;
	return true;
}

bool
operator>(const vec3 &vecA, const vec3 &vecB){ 
	if (vecA[0] < vecB[0]) return false;
	if (vecA[1] < vecB[1]) return false;
	if (vecA[2] < vecB[2]) return false;
	return true;
}

float AABB::intersect(Ray& ray){
	
	if (ray.origin < aabbmax && ray.origin > aabbmin){
		return EPSILON; // always first if inside	
	}
	
	
	float far = FLT_MAX;
	float near = FLT_MIN;
	
	if (!intersect1D(ray.origin[0],ray.direction[0],aabbmin[0],aabbmax[0],near,far)) return false;
	if (!intersect1D(ray.origin[1],ray.direction[1],aabbmin[1],aabbmax[1],near,far)) return false;
	if (!intersect1D(ray.origin[2],ray.direction[2],aabbmin[2],aabbmax[2],near,far)) return false;
	
	return near+EPSILON;
}






