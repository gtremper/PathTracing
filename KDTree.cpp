/* Implementation of KDTree acceleration structure */

#include <iostream>
#include <vector>

#include "Shapes.h"
#include "Intersection.h"
#include "KDTree.h"

using namespace std;

#define TREELIMIT 5
#define EPSILON 0.000001

/***  TREENODE  ***/
struct ShapeSorter {
	ShapeSorter(int ax){this->axis = ax;}
	bool operator() (Shape* i, Shape* j) {
		return i->aabb.center[axis] < j->aabb.center[axis];
	}
	int axis;
};

TreeNode::TreeNode(vector<Shape*>& prims, AABB& bigBox){
	this->aabb = bigBox;
	
	vec3 diff = bigBox.max - bigBox.min;
	if (diff[0]>diff[1] && diff[0]>diff[1]){
		this->axis = 0;
	} else if (diff[1]>diff[2]){
		this->axis = 1;
	} else {
		this->axis = 2;
	}
	
	ShapeSorter s(axis);
	sort(prims.begin(), prims.end(),s);
	double split = prims[prims.size()/2]->aabb.center[axis];
	
	vector<Shape*> leftPrims;
	vector<Shape*>::iterator it=prims.begin();
	for(;it!=prims.end(); ++it){
		if((*it)->aabb.min[axis]<(split+EPSILON)){
			leftPrims.push_back(*it);
		}
	}
	
	vector<Shape*> rightPrims;
	for(it=prims.begin(); it!=prims.end(); ++it){
		if((*it)->aabb.max[axis]>(split-EPSILON)){
			rightPrims.push_back(*it);
		}
	}
	
	unsigned int threshold = (prims.size()/TREELIMIT) * (2*TREELIMIT-1);
	unsigned int rsize = rightPrims.size();
	unsigned int lsize = leftPrims.size();
	bool above_tresh = (lsize + rsize) > threshold;
	
	if( above_tresh || prims.size()==lsize || prims.size()==rsize || prims.size()<=TREELIMIT){
		left = NULL;
		right = NULL;
		primatives = prims;
	}else{
		AABB leftAABB = AABB(aabb);
		leftAABB.max[axis] = split;
		
		AABB rightAABB = AABB(aabb);
		rightAABB.min[axis] = split;
		
		left = new TreeNode(leftPrims, leftAABB);
		right = new TreeNode(rightPrims, rightAABB);
	}
	
}

TreeNode::~TreeNode(){
	if(left){
		delete left;
		delete right;
	}
}

Intersection TreeNode::intersect(Ray& ray){
	if(!left){
		return Intersection(primatives,ray);
	}
	
	double hitLeft = left->aabb.intersect(ray,this->axis);
	double hitRight = right->aabb.intersect(ray,this->axis);
	
	if(!hitLeft && !hitRight) return Intersection();
	
	if (!hitLeft) {
		return right->intersect(ray);
	}
	if (!hitRight) {
		return left->intersect(ray);
	}
	
	Intersection hit;
	if (hitLeft > hitRight){
		hit = right->intersect(ray);
		if (!hit.primative){
			hit = left->intersect(ray);
		}
		return hit;
	} else {
		hit = left->intersect(ray);
		if (!hit.primative){
			hit = right->intersect(ray);
		}
		return hit;
	}
		
}
