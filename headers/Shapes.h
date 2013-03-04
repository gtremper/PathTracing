/* Class definition for a shape */
#ifndef SHAPE_H
#define SHAPE_H

#include "Transform.h"
#include "Ray.h"

class AABB {
public:
	AABB(){};
	AABB(vec3&,vec3&);
	vec3 aabbmax;
	vec3 aabbmin;
	float intersect(Ray&);
	vec3 center;
};

class Shape {
	public:
		virtual float intersect(Ray&)=0;
		virtual vec3 getNormal(vec3&)=0;
	//	virtual vec3 getTexture(vec3&){};
		
		AABB aabb;
		bool hasTexture;
		
		/* material properties */
		vec3 ambient;
		vec3 diffuse;
		vec3 specular;
		float shininess;
		vec3 emission;
		float indexofrefraction;
		float refractivity;
};

class Sphere : public Shape {
	public:
		Sphere(mat4);
		float intersect(Ray&);
		vec3 getNormal(vec3&);
	
		mat4 mv;
		mat4 inv;	
};

class Triangle : public Shape {
	public:
		Triangle(){};
		Triangle(vec3,vec3,vec3);
		float intersect(Ray&);
		virtual vec3 getNormal(vec3&);
	//	vec3 getTexture(vec3&);
	
		vec3 p0;
		vec3 p1;
		vec3 p2;
		vec3 n0;
};

class NormTriangle : public Triangle {
public:
	NormTriangle(vec3,vec3,vec3,vec3,vec3,vec3);
	vec3 getNormal(vec3&);

	vec3 n1;
	vec3 n2;	
};

#endif 