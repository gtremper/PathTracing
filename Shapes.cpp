/* Class definitions for Triangles and Spheres */
#include <iostream>
#include <vector>

#include "Shapes.h"
#include "Intersection.h"
#include "KDTree.h"

#define EPSILON 0.000001

using namespace std;

/***  RAY  ***/
Ray::Ray(const vec3& o,const vec3& d){
	origin = o;
	direction = glm::normalize(d);
}

vec3 Ray::getPoint(double t) {
	return origin + t*direction;
}

/***  INTERSECTION  ***/
Intersection::Intersection(vector<Shape*>& objects, Ray& ray) {
	double min_t = DBL_MAX;
	primative = NULL;
	for(vector<Shape*>::iterator prim=objects.begin(); prim!=objects.end(); prim++) {
		double t = (*prim)->intersect(ray);
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
	double minx = min(min(p0[0],p1[0]),p2[0]);
	double maxx = max(max(p0[0],p1[0]),p2[0]);
	double miny = min(min(p0[1],p1[1]),p2[1]);
	double maxy = max(max(p0[1],p1[1]),p2[1]);
	double minz = min(min(p0[2],p1[2]),p2[2]);
	double maxz = max(max(p0[2],p1[2]),p2[2]);
	vec3 minvec = vec3(minx,miny,minz);
	vec3 maxvec = vec3(maxx,maxy,maxz);
	aabb = AABB(minvec, maxvec);
}


double Triangle::intersect(Ray& ray){
	vec3 col1 = p1-p0;
	vec3 col2 = p2-p0;
	mat3 M = mat3(col1, col2, -ray.direction);
	double det = glm::determinant(M);
	if (det<EPSILON && det>-EPSILON) return -1.0;
	
	M[0] = ray.origin-p0;
	double alpha = glm::determinant(M)/det;
	if (alpha<0.0 || alpha>1.0) return -1.0;
	M[0] = col1;
	M[1] = ray.origin-p0;

	double beta = glm::determinant(M)/det;
	if (beta<0.0 || beta+alpha>1.0) return -1.0;
	M[1] = col2;
	M[2] = ray.origin-p0;
	return glm::determinant(M)/det;
}

vec3 Triangle::getNormal(vec3& hit){
	return n0;
}


double Triangle::getSubtendedAngle(const vec3& origin) {
	vec3 a = p0 - origin;
	vec3 b = p1 - origin;
	vec3 c = p2 - origin;
	
	double det = glm::determinant( mat3(a,b,c) );
	
	double a1 = glm::length(a);
	double b1 = glm::length(b);
	double c1 = glm::length(c);
	
	double result = a1*b1*c1 + glm::dot(a,b)*c1 + 
				glm::dot(a,c)*b1 + glm::dot(b,c)*a1;
					
	result = atan2(det, result);
	if (result < 0.0) result += M_PI;
	result *= 2.0;
	return result;
}

/*
hit is the point hit on the triangle
light is the light you want to shade with
tree is the KDtree
num_samples is the root of the number
*/
vec3 Triangle::shade(Intersection& hit, Shape* light, 
						TreeNode* tree,double num_samples){
	double subdivide = 1.0 / num_samples;
	double numHits = 0.0;
	vec3 thisNormal = getNormal(hit.point);
	for(double a=0; a<num_samples; a+=1) {
		for(double b=0; b<num_samples; b+=1) {
			//calculate statified random point on triangle
			double rand1 = ((double)rand()/(double)RAND_MAX) * subdivide;
			double rand2 = ((double)rand()/(double)RAND_MAX) * subdivide;
			vec3 dir = p0;
			dir += (a*subdivide + rand1) * (p1-p0);
			dir += (b*subdivide + rand2) * (p2-p0);
			dir = glm::normalize(dir - hit.point);
			Ray ray = Ray(hit.point + EPSILON * thisNormal, dir);
			Intersection light_hit = tree->intersect(ray);
			if (light_hit.primative == this) {
				numHits += 1.0;
			}
		}
	}
	
	vec3 lightNormal = light->getNormal(aabb.center);
	
	vec3 direction = glm::normalize(this->aabb.center - hit.point);
	vec3 shade = max(0.0,glm::dot(thisNormal,direction)) * diffuse;

	vec3 half = glm::normalize(hit.sourceDirection+direction);
	double phong = pow( max(0.0,glm::dot(half,thisNormal)) , shininess);
	shade += phong * specular;
	shade *= light->emission;
	
	double dist = glm::distance(hit.point, this->aabb.center);
	dist *= dist; // square distance
	double cos_weight = glm::dot(thisNormal,direction);
	cos_weight *= glm::dot(lightNormal, -direction);
	cos_weight /= dist;
	
	double frac = numHits * subdivide*subdivide; //fraction shadow rays hit
	return frac * cos_weight * shade;
}

/*
vec3 Triangle::getTexture(vec3& hit){
	mat2 M = mat2(p1[0]-p0[0], p1[1]-p0[1], p2[0]-p0[0], p2[1]-p0[1]);
	double det = glm::determinant(M);
	M[0][0] = hit[0]-p0[0];
	M[0][1] = hit[1]-p0[1];
	double alpha = glm::determinant(M)/det;//p0 to p1
	M[0][0] = p1[0]-p0[0];
	M[0][1] = p1[1]-p0[1];
	M[1][0] = hit[0]-p0[0];
	M[1][1] = hit[1]-p0[1];
	double beta = glm::determinant(M)/det;//p0 to p2
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
	double det = glm::determinant(M);
	M[0][0] = hit[0]-p0[0];
	M[0][1] = hit[1]-p0[1];
	double beta = glm::determinant(M)/det;
	M[0][0] = p1[0]-p0[0];
	M[0][1] = p1[1]-p0[1];
	M[1][0] = hit[0]-p0[0];
	M[1][1] = hit[1]-p0[1];
	double gamma = glm::determinant(M)/det;
	return glm::normalize((1-beta-gamma)*n0 + beta*n1 + gamma*n2);
}

/***  SPHERE  ***/
Sphere::Sphere(mat4 trans){
	mv = trans;
	inv = glm::inverse(mv);
	
	/*Bounding Box*/
	mat4 S = mat4(1.0);
	S[3][3] = -1.0;
	mat4 R = mv*S*glm::transpose(mv);
	double minx = (R[0][3]+sqrt(R[0][3]*R[0][3] -R[3][3]*R[0][0]))/R[3][3];
	double maxx = (R[0][3]-sqrt(R[0][3]*R[0][3] -R[3][3]*R[0][0]))/R[3][3];
	double miny = (R[1][3]+sqrt(R[1][3]*R[1][3] -R[3][3]*R[1][1]))/R[3][3];
	double maxy = (R[1][3]-sqrt(R[1][3]*R[1][3] -R[3][3]*R[1][1]))/R[3][3];
	double minz = (R[2][3]+sqrt(R[2][3]*R[2][3] -R[3][3]*R[2][2]))/R[3][3];
	double maxz = (R[2][3]-sqrt(R[2][3]*R[2][3] -R[3][3]*R[2][2]))/R[3][3];
	
	vec3 minvec = vec3(minx,miny,minz);
	vec3 maxvec = vec3(maxx,maxy,maxz);
	aabb = AABB(minvec, maxvec);
}



double Sphere::intersect(Ray& ray) {
	vec3 direction = glm::normalize(vec3(inv * vec4(ray.direction,0.0)));
	vec3 origin =vec3(inv * vec4(ray.origin,1.0));
	
	double b = 2.0f * glm::dot(direction, origin);
	double c = glm::dot(origin,origin) - 1.0;
	double det = b*b - 4.0f*c;
	if (det<0.0) return -1.0;
	det = sqrt(det);
	
	double t1 = (-b+det)*0.5;
	double t2 = (-b-det)*0.5;
	
	if (t1<0.0 && t2<0.0) return -1.0;
	
	if (t2>0.0f) {
		vec4 hit = mv * vec4(origin+t2*direction,1.0);
		return glm::distance(ray.origin,vec3(hit));
	} else { //t1 is closer
		vec4 hit = mv * vec4(origin+t1*direction,1.0);
		return glm::distance(ray.origin,vec3(hit));
	}
}

vec3 Sphere::getNormal(vec3& hit){
	return glm::normalize(vec3(glm::transpose(inv)*inv*vec4(hit,1.0)));
}

double Sphere::getSubtendedAngle(const vec3& origin) {
	cout << "THIS IS NOT IMPLEMENTED" << endl;
	exit(1);
}

vec3 Sphere::shade(Intersection& hit, Shape* light, 
						TreeNode* tree,double num_samples){
	cout << "THIS IS NOT IMPLEMENTED" << endl;
	exit(1);
}


/***  AABB  ***/
AABB::AABB(vec3& minarg, vec3& maxarg) {
	this->aabbmin = minarg;
	this->aabbmax = maxarg;
	this->center = (minarg+maxarg)/2.0;
	
}

inline
bool intersect1D(double start, double dir, double axisMin, double axisMax, double& near, double& far){
	// Parallel
	if(dir<EPSILON && dir>-EPSILON){
		return (start>axisMin) && (start<axisMax);
	}
	
	//intersection parameters
	double t0 = (axisMin-start)/dir;
	double t1 = (axisMax-start)/dir;
	
	if(t0>t1){
		double temp = t1;
		t1 = t0;
		t0 = temp;
	}
	
	near = max(t0,near);
	far = min(t1,far);
	
	if(near>far || far<0.0) return false;
	
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

double AABB::intersect(Ray& ray){
	
	if (ray.origin < aabbmax && ray.origin > aabbmin){
		return EPSILON; // always first if inside	
	}
	
	
	double far = DBL_MAX;
	double near = DBL_MIN;
	
	if (!intersect1D(ray.origin[0],ray.direction[0],aabbmin[0],aabbmax[0],near,far)) return false;
	if (!intersect1D(ray.origin[1],ray.direction[1],aabbmin[1],aabbmax[1],near,far)) return false;
	if (!intersect1D(ray.origin[2],ray.direction[2],aabbmin[2],aabbmax[2],near,far)) return false;
	
	return near+EPSILON;
}






