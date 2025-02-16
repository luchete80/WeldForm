//#include "Mesh.h"
#include <iostream>
using namespace std;

namespace SPH {

	
TriMesh::TriMesh(){
	
  m_v = 0.;
  m_w = 0.;
  dimension = 3;
	
}

//TODO: CHANGE TRIMESH NAME
TriMesh::TriMesh(NastranReader &nr, bool flipnormals){
  dimension = nr.dim;
  //Insert nodes
  for (int n=0;n<nr.node_count;n++){
    //if (!flipnormals)
      node.Push(new Vec3_t(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]));
    // else 
      // node.Push(new Vec3_t(nr.node[3*n+1],nr.node[3*n],nr.node[3*n+2]));
    
		node_v.Push(new Vec3_t(0.,0.,0.));
  }
  cout << "Generated "<<node.Size()<< " trimesh nodes. "<<endl;
  //cout << "Normals"<<endl;
  cout << "Writing elements..."<<endl;
  for (int e=0;e<nr.elem_count;e++){
    double cz = nr.elcon[3*e+2];
    if (dimension == 2) cz = 0;
    if (!flipnormals)   element.Push(new Element(nr.elcon[3*e],nr.elcon[3*e+1],cz));		  
    else                element.Push(new Element(nr.elcon[3*e+1],nr.elcon[3*e],cz));		        
    Vec3_t v;
		if (dimension ==3) v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]] + *node[nr.elcon[3*e+2]] ) / 3. ;
    else               v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]])  / 2. ;
    element[e] -> centroid = v;
    //TODO: CHANGE FOR CALCNORMALS
    if (dimension==3){
      Vec3_t v1, v2;
      //In COUNTERCLOCKWISE
      v1 = *node[nr.elcon[3*e+1]] - *node[nr.elcon[3*e]];
      v2 = *node[nr.elcon[3*e+2]] - *node[nr.elcon[3*e]];
      element[e] ->normal = cross (v1,v2);

      element[e] ->normal /= Norm(element[e] ->normal);
      //cout << "v1 "<< v1<< ", v2 " <<v2<< ", normal "<<element[e]->normal <<endl;
    } else { //See calc normals
        Vec3_t u = *node [element[e]->node[1]] - *node [element[e]->node[0]];
        v[0] = -u[1];
        v[1] =  u[0];
        v[2] =  0.0;
        element[e] -> normal = v/norm(v);
    }
  }
  cout << "Generated "<<element.Size()<< " trimesh elements. "<<endl;  
  
  m_v = 0.;
  m_w = 0.;
}

Element::Element(const int &n1, const int &n2, const int &n3){
	
	//centroid = Vec3_t();
	node[0] = n1; node [1] = n2; node[2] = n3;
	
}

void TriMesh::CalcCentroids(){
	
	for (int e=0;e<element.Size();e++){
		element[e]-> centroid = 0.;
    for (int i=0;i<dimension;i++)
      element[e]-> centroid += *node[element[e]->node[i]];
    element[e]-> centroid/= dimension; 
	}
}

// TODO: extend to all dirs
inline void TriMesh::AxisPlaneMesh(const int &axis, bool positaxisorent, const Vec3_t p1, const Vec3_t p2,  const int &dens){
	int elemcount = dens * dens;
  
  if (dimension == 2) elemcount = dens; 
	
	double x1,x2,x3;
	double l1,l2;
	Vec3_t p = p2-p1;
	int dir[3];
	if 			(axis == 0 )	{dir[0] = 1; dir[1] = 2;}
	else if (axis == 1 )	{dir[0] = 0; dir[1] = 2;}
	else									{dir[0] = 0; dir[1] = 1;}
	
	dir [2] = axis; //dir2 is which remains constant
	
	x3 = p1(dir[2]);

	x2=p1(dir[1]); 
	double dl = p(dir[0])/dens;	//Could be allowed 2 diff densities
	//cout <<"dens: "<<dens<<endl;
	//Plane is in 0 and 1 dirs
	int test =dens+1;
  
  if (dimension == 2) test = 1;
  
	for (int j=0; j<test; j++) {
		x1 = p1(dir[0]);
		for (int i=0; i<dens+1; i++){
			Vec3_t v;
			v(dir[0])=x1;v(dir[1])=x2;v(dir[2])=x3;
			//cout << "i,j" << i << ", " << j<<endl; 
			//node.Push(new Vec3_t(x1,x2,x3));
			node.Push(new Vec3_t(v(0),v(1),v(2)));
			node_v.Push(new Vec3_t(0.,0.,0.));
			//cout << "x y : "<<v(0) << ", "<<v(1)<<", "<<v(2)<<endl;
			x1+=dl;
		}
		x2+=dl;
	}
  cout << "Created "<<node.Size()<< " nodes "<<endl;

	int n[4];
	int el =0;
	int i;
  cout << "Creating elements "<<endl;
  if (dimension == 3){
    for (size_t j = 0 ;j  < dens; j++ ) {
          // cout <<"j, dens" <<j<<", "<<dens<<endl;
          // cout <<"j<dens"<< (j  < dens)<<endl;
      for ( i = 0; i < dens; i++ ){
          // cout <<"i, dens" <<i<<", "<<dens<<endl;
          // cout <<"i <dens"<< (i  < dens)<<endl;
          n[0] = (dens + 1)* j + i; 		n[1] = n[0] + 1; 
          n[2] = (dens + 1)* (j+1) + i; n[3] = n[2] + 1;
        //cout <<" jj" << jj<<endl;
        int elcon[2][3];	// TODO: check x, y and z normals and node direction 
                          // For all plane orientations
        //If connectivity  is anticlockwise normal is outwards
        if (positaxisorent) {
          elcon[0][0] = n[0];elcon[0][1] = n[1];elcon[0][2] = n[2];
          elcon[1][0] = n[1];elcon[1][1] = n[3];elcon[1][2] = n[2];
        } else {
          elcon[0][0] = n[0];elcon[0][1] = n[2];elcon[0][2] = n[1];
          elcon[1][0] = n[1];elcon[1][1] = n[2];elcon[1][2] = n[3];				
        }
        //cout << "elnodes"<<endl;
        for ( int e= 0; e<2;e++) { // 2 triangles
          element.Push(new Element(elcon[e][0],elcon[e][1],elcon[e][2]));		
          //cout << "Element "<< el <<": ";
          // for (int en = 0 ; en<3; en++) cout << elcon[e][en]<<", ";
          // cout <<endl;
          
          Vec3_t v = ( *node[elcon[e][0]] + *node[elcon[e][1]] + *node[elcon[e][2]] ) / 3. ;
          element[el] -> centroid = v; 
          //cout << "Centroid" << element[el] -> centroid << endl;
          el++;
        }
      }// i for
      
    }
  } else {
    for ( i = 0; i < dens; i++ ){
          n[0] = i; 		n[1] = n[0] + 1; 
        //cout <<" jj" << jj<<endl;
        int elcon[2];	// TODO: check x, y and z normals and node direction 
                          // For all plane orientations
        if (positaxisorent) {
          elcon[0] = n[0];elcon[1] = n[1];
        } else {
          elcon[0] = n[1];elcon[1] = n[0];		
        }

        element.Push(new Element(elcon[0],elcon[1],0));		
        //cout << "Element "<< el <<": ";
        // for (int en = 0 ; en<3; en++) cout << elcon[e][en]<<", ";
        // cout <<endl;
        
        Vec3_t v = ( *node[elcon[0]] + *node[elcon[1]] ) / 2. ;
        element[el] -> centroid = v; 
        //cout << "Centroid" << element[el] -> centroid << endl;
        el++;                
    }  
  }
	///////////////////////////////////////////
	//// MESH GENERATION END
	cout << "Creating normals"<<endl;
  CalcNormals();
  if (positaxisorent){
    for (int e = 0; e < element.Size(); e++)
      element[e] -> normal= -element[e] -> normal;
  }  
	// for (int e = 0; e < element.Size(); e++){ 
		// double f=-1.;
		// if (positaxisorent) f= 1.;
		// element[e] -> normal (axis) = f;
	// }

  cout << "Created Mesh with "<< node.size()<< " nodes. "<<endl;
  if (node.size() == 0)
    throw new Fatal("ATTENTION! Check mesh generation");
}

//This is done once, Since mesh is rigid
//Calculate radius and plane coefficient
inline void TriMesh::CalcSpheres(){
	double max;
	for (int e = 0; e < element.Size(); e++){ 
		max = 0.;
		Vec3_t rv;
		for (int n = 0 ;n < dimension; n++){
			rv = *node [element[e]->node[n]] - element[e] -> centroid;
			if (norm(rv) > max) max = norm(rv);
			element[e]-> nfar = n;
		}
		element[e]-> radius = max;	//Fraser Eq 3-136
	}
	UpdatePlaneCoeff();
	
}

inline void TriMesh::UpdatePlaneCoeff(){
	//Update pplane
	for (int e = 0; e < element.Size(); e++) 
		element[e]-> pplane = dot(*node [element[e] -> node[element[e] ->nfar]],element[e] -> normal);

}
// TODO: PARALELIZE
inline void TriMesh::CalcNormals(){
	Vec3_t u, v, w;
  
  if (dimension==3){
    for (int e = 0; e < element.Size(); e++) {
        u = *node [element[e]->node[1]] - *node [element[e]->node[0]];
        v = *node [element[e]->node[2]] - *node [element[e]->node[0]];
        w = cross(u,v);
        element[e] -> normal = w/norm(w);
        //Fraser Eqn 3.34
        //Uj x Vj / |UjxVj|
    }
  } else {//ROTATE COUNTERCLOCKWISE (SURFACE BOUNDARY IS SORROUNDED CLOCKWISE to outer normal)
    /////// i.e. : x.positive line has y positive normal
      for (int e = 0; e < element.Size(); e++) {
        u = *node [element[e]->node[1]] - *node [element[e]->node[0]];
        v[0] = -u[1];
        v[1] =  u[0];
        v[2] =  0.0;
        element[e] -> normal = v/norm(v);
        //cout << "Element " << e << " normal "<<element[e] -> normal<<endl;
      }
      //x2=cosβx1−sinβy1
      //y2=sinβx1+cosβy1
      //Sin (PI/2.) = -1
      //Cos (PI/2.) = 0,
  }
}

// ATENTION: THIS OVERRIDES AUTO UPDATE WITH VELOCITIES 
// IF RIGID TRANSLATION AND ROTATION IS USED; THIS IS NOT PREFERRED
inline void TriMesh::ApplyConstVel(const Vec3_t &v){
		for (int n=0;n<node.Size();n++)
			*node_v[n] = v;
}

inline void TriMesh::CalcCentroidVelFromNodes(){
	
	
}

//UPDATES MESH NODES POSITIONS, VELOCITIES, NORMALS AND PPLANE
//ALL FROM RIGID TRANSLATION AND ROTATION
//THIS USES m_v and m_w members
inline void TriMesh::Update(const double &dt){
	//Seems to be More accurate to do this by node vel
	//This is used by normals
  Vec3_t min = 1000.;
  Vec3_t max = -1000.;
	for (int n=0;n<node.Size();n++){
    Vec3_t vr 	= cross(m_w, *node[n]);
    *node_v[n] = m_v + vr;
    for (int i=0;i<dimension;i++) {
      if      ((*node[n])(i) < min(i)) min[i] = (*node[n])(i);
      else if ((*node[n])(i) > max(i)) max[i] = (*node[n])(i);
    } 
		*node[n] += (*node_v[n])*dt;
	}
  
  //cout << "Min Max Node pos" << min<< "; " <<max<<endl;
  
  CalcCentroids();
  CalcNormals();        //From node positions
  UpdatePlaneCoeff();   //pplane
}
 
//To use in contact intersection
inline void TriMesh::Move(const Vec3_t &v){
	//Seems to be More accurate to do this by node vel
	//This is used by normals

	for (int n=0;n<node.Size();n++){
		*node[n] += v;
	} 

  CalcCentroids();
  CalcNormals();        //From node positions
  UpdatePlaneCoeff();   //pplane
  
}

inline void TriMesh::Scale(const double &f){
	//Seems to be More accurate to do this by node vel
	//This is used by normals

	for (int n=0;n<node.Size();n++){
		*node[n] *= f;
	} 
  cout << "calc centroids"<<endl;
  CalcCentroids();
  cout << "calc normals"<<endl;
  CalcNormals();        //From node positions
  cout << "generate plane coeffs"<<endl;
  //UpdatePlaneCoeff();   //pplane
}
 
  
  
};
