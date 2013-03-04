#ifndef SCENE_H  
#define SCENE_H

#include "KDTree.h"

class Scene {
	public:
		vec3 w;
		vec3 u;
		vec3 v;
		int width;
		int height;
		int maxdepth;
		std::string filename;
		float fovy;
		float fovx;
		vec3 eye;
		std::vector<Shape*> objects;
		AABB sceneAABB;
		TreeNode* KDTree;
		
		vec3 ambient;
		vec3 diffuse;
		vec3 specular;
		float shininess;
		vec3 emission;
		float indexofrefraction;
		float refractivity;
		
		int antialias;
		int shadowrays;
		float lightradius;
		
		std::vector<Light*> lights;
		float constant;
		float linear;
		float quadratic;
		
		Scene(char*);
		~Scene();
		Ray castEyeRay(float,float);
		void setCoordinateFrame(vec3&,vec3&);
		void parseLine(std::string,std::stack<mat4>&,
			std::vector<vec3>&,std::vector<vec3>&,std::vector<vec3>&);
		void parse(char*);
		void updateAABB(vec3&);
};


#endif