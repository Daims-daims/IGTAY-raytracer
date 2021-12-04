
#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <stdio.h>

#include <glm/gtc/epsilon.hpp>

/// acne_eps is a small constant used to prevent acne when computing
/// intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {
  //! \todo : compute intersection of the ray and the plane object
  float denom = dot(ray->dir,obj->geom.plane.normal);
  if(denom==0){
    return false;
  }
  float num = -(dot(ray->orig,obj->geom.plane.normal)+obj->geom.plane.dist);
  float res = num/denom;
  if(ray->tmax< res || ray->tmin > res ){
    return false;
  }
  intersection->position = ray->orig + (res) * ray->dir;
  intersection->normal = normalize(obj->geom.plane.normal);
  intersection->mat=&obj->mat;
  return true;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {

  //! \todo : compute intersection of the ray and the sphere object
  vec3 O_C =  ray->orig - obj->geom.sphere.center;
  float R_2 = obj->geom.sphere.radius;
  float a = 1;
  float b = 2*dot(ray->dir,O_C);
  float c = dot(O_C,O_C)-R_2*R_2;
  float discr = b*b - 4 * a * c;
  float res;
  if ( discr < 0 ){
    return false;
  }
  if ( discr > 0 ){
    float x_1 = (-b - sqrt(discr))/(2*a);
    float x_2 = (-b + sqrt(discr))/(2*a);
    float res_tmp_1 = min(x_1,x_2);
    float res_tmp_2 = max(x_1,x_2);
    if(ray->tmax> res_tmp_1 && ray->tmin < res_tmp_1){
      res=res_tmp_1;
    }
    else if(ray->tmax> res_tmp_2 && ray->tmin < res_tmp_2){
      res=res_tmp_2;
    }
    else{
      return false;
    }
  }
  else{
    float x = (-b /(2*a));
    if(ray->tmax< x || ray->tmin > x){
      return false;
    }
    res=x;
  }
  intersection->position=ray->orig+ray->dir*res;
  intersection->normal=normalize(intersection->position - obj->geom.sphere.center );
  intersection->mat=&obj->mat;
  return true;
}

bool x1_leq_x2(Intersection *intersection,Intersection *intersectionbis,Ray* ray){
  float d1 = dot(intersection->position-ray->orig,intersection->position-ray->orig);
  float d2 = dot(intersectionbis->position-ray->orig,intersectionbis->position-ray->orig);
  return d1<=d2;
}

void echangerIntersection(Intersection* nvIntersection, Intersection* oldIntersection){
  nvIntersection->normal=oldIntersection->normal;
  nvIntersection->position=oldIntersection->position;
  nvIntersection->mat=oldIntersection->mat;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  bool hasIntersection = false;
  bool tmpIntersection = false;
  int objectCount = scene->objects.size();
  Object* obj;
  Intersection intTmp;


//!\todo loop on each object of the scene to compute intersection
// and find closer object from point of view 
  for(int i = 0 ; i< objectCount; i++){
    obj = scene->objects[i];
    if(obj->geom.type==SPHERE && intersectSphere(ray,&intTmp,obj)) tmpIntersection = true ;
    if(obj->geom.type==PLANE && intersectPlane(ray,&intTmp,obj)) tmpIntersection = true ;
    if(tmpIntersection){
      if(!hasIntersection) echangerIntersection(intersection,&intTmp);
      else if(x1_leq_x2(&intTmp,intersection,ray)) echangerIntersection(intersection,&intTmp);
      hasIntersection=true;
    }
    tmpIntersection=false;
  }
  return hasIntersection;
}

bool intersectSceneOmbre(const Scene *scene, Ray *ray){
  int objectCount=scene->objects.size();
  Intersection intTmp;
  Object* obj;
  for(int i = 0 ; i< objectCount; i++){
    obj = scene->objects[i];
    if(obj->geom.type==SPHERE && intersectSphere(ray,&intTmp,obj)) return true ;
    if(obj->geom.type==PLANE && intersectPlane(ray,&intTmp,obj)) return true ;
  }
  return false;
}

/* ---------------------------------------------------------------------------
 */
/*
 *	The following functions are coded from Cook-Torrance bsdf model
 *description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba
 *renderer)
 */

// Shadowing and masking function. Linked with the NDF. Here, Smith function,
// suitable for Beckmann NDF
float RDM_chiplus(float c) { return (c > 0.f) ? 1.f : 0.f; }

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */

float RDM_Beckmann(float NdotH, float alpha) {
  ;
  float tan = (1-(NdotH*NdotH))/(NdotH*NdotH);
  float denom = M_PI * alpha*alpha * NdotH*NdotH*NdotH*NdotH;
  float num = exp(-(tan)/(alpha*alpha));
  return RDM_chiplus(NdotH)*num/denom;
}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {
  float sin_2 = (extIOR/intIOR)*(extIOR/intIOR)*(1-(LdotH*LdotH));
  if(sin_2 > 1 ){ return 1.f ;};
  float cos_t = sqrt(1-sin_2);
  float cos_i = LdotH; 
  float R_S_h = (extIOR * cos_i - intIOR * cos_t);
  float R_S_b = (extIOR * cos_i + intIOR * cos_t);
  float R_P_h = (extIOR * cos_t - intIOR * cos_i);
  float R_P_b = (extIOR * cos_t + intIOR * cos_i);
  float R_s = (R_S_h*R_S_h ) / (R_S_b*R_S_b);
  float R_p = (R_P_h*R_P_h) / (R_P_b*R_P_b);
  float addition = R_s+R_p;
  return addition/2.f;  

}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {
float tan = sqrt(1-(DdotN*DdotN))/DdotN;
float k = DdotH/DdotN;
float b = 1/(alpha*tan);
if( b>=1.6) { return RDM_chiplus(k);}
float num = 3.535*b + 2.181*b*b;
float deno = 1 + 2.276*b + 2.577*b*b;
return RDM_chiplus(k)*num/(deno);

}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                float alpha) {
  return (RDM_G1(LdotH,LdotN,alpha)*RDM_G1(VdotH,VdotN,alpha));
}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, Material *m) {
  color3 spec = m-> specularColor;
  float alpha = m->roughness;
  float D = RDM_Beckmann(NdotH,alpha);
  float F = RDM_Fresnel(LdotH,1.f,m->IOR);
  float G = RDM_Smith(LdotH,LdotN,VdotH,VdotN,alpha);
  float num = D *F *G;
  float deno = 4 * LdotN * VdotN;
  float mult = num/deno;
  return spec*mult;

}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {

  return m->diffuseColor/(float)M_PI;

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                Material *m) {

  return (RDM_bsdf_d(m)+RDM_bsdf_s(LdotH,NdotH,VdotH,LdotN,VdotN,m));

}
/*
float RDM_Beckmann(float NdotH, float alpha) {


  //! \todo compute Beckmann normal distribution
  if(NdotH<=0){
    return 0.f;
  }
  float alpha_2 = alpha * alpha;
  float cos_2 = NdotH * NdotH;
  float tan_2 = (1-cos_2)/(cos_2);
  float nume = exp(-tan_2/alpha_2);
  float deno = M_PI * alpha_2 * NdotH * NdotH * NdotH * NdotH;
  //if (nume/deno != 0 ) printf("%f / %f = %f \n",-tan_2/alpha_2,deno,nume/deno);
  return nume/deno;

}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {

  //! \todo compute Fresnel term
  // ext = n1 int = n2

  float cos_2 = LdotH*LdotH;
  float sin_2 = (extIOR/intIOR)*(extIOR/intIOR)*(1-cos_2); 
  if(sin_2>1) return 1.f;
  float cost = sqrt(1-sin_2);
  float R_s_h = (extIOR * LdotH - intIOR * cost);
  float R_s_b = (extIOR * LdotH + intIOR * cost);
  float R_p_h = (extIOR * cost - intIOR * LdotH);
  float R_p_b = (extIOR * cost + intIOR * LdotH);
  float R_s = (R_s_h*R_s_h)/(R_s_b*R_s_b);
  float R_p = (R_p_h*R_p_h)/(R_p_b*R_p_b);
  return (R_s+R_p)/2;
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {

  //! \todo compute G1 term of the Smith fonction
  float tan= (sqrt(1-(DdotN*DdotN)))/DdotN;
  float b = 1/(alpha*tan);
  float k = DdotH/DdotN;
  float lambda=0;
  //if (b<1)printf(" b : 1 / %f * %f = %f \n",alpha,tan,b);
  if(k>0) lambda = 1;
  else return 0;
  if ( b >1.6) return lambda;
  float h = 3.535*b + 2.181*b*b;
  float b2 = 1 + 2.276*b + 2.577*b*b;
  //printf("G1 : %f/%f=%f\n",h,b2,h/b2);
  return h/b2;

}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                float alpha) {
  //! \todo the Smith fonction
  //if (VdotH !=0 ) printf(" Smith : L.H = %f - L.N = %f - V.H = %f - V.N = %f \n",LdotH,LdotN,VdotH,VdotN);
  float V = RDM_G1(VdotH,VdotN,alpha);
  float L = RDM_G1(LdotH,LdotN,alpha);
  //printf("G : %f * %f = %f \n",V,L,V*L);
  return L*V;

}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, Material *m) {

  //! \todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G
  //! = RDM_Smith
  float alpha = m->roughness;
  float deno = 4*VdotN*LdotN;
  float D = RDM_Beckmann(NdotH,alpha);
  float F = RDM_Fresnel(LdotH,1.f,m->IOR);
  float G = RDM_Smith(LdotH,LdotN,VdotH,VdotN,alpha);
  float num = D*F*G;
  float mult = num*deno;
  //if(D>.0001) printf(" res : %f-%f-%f = %f\n",D,F,G,mult);
  return m->specularColor * mult;

}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {

  //! \todo compute diffuse component of the bsdf
  return m->diffuseColor*(1.f/(float) M_PI);

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                Material *m) {

  //! \todo compute bsdf diffuse and specular term
  color3 d = RDM_bsdf_d(m);
  color3 s = RDM_bsdf_s(LdotH,NdotH,VdotH,LdotN,VdotN,m);
  color3 res = d+s;
  //printf("d : %f - %f - %f \n s : %f - %f - %f\n d + s : %f - %f - %f\n", d[0],d[1],d[2],s[0],s[1],s[2],res[0], res[1], res[2]);
  return (d+s);

}*/


color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat) {
  color3 ret = color3(0.f);

//! \todo compute bsdf, return the shaded color taking into account the
//! lightcolor
  vec3 h_h = normalize(l+v);
  if((float)dot(l,n)<=0){
    return ret;
  }
  color3 BSDF= RDM_bsdf(dot(l,h_h),dot(n,h_h),dot(v,h_h),dot(l,n),dot(v,n),mat);
  ret = (float)dot(l,n) * lc * BSDF ; 
  return ret;
}

//! if tree is not null, use intersectKdTree to compute the intersection instead
//! of intersect scene
/*
color3 calculShade(Scene* scene,Intersection* intersection, vec3 v , Lights l_s){
  int ln = l_s.size(); 
  color3 ret = color3(0,0,0);
  vec3 l= vec3 ( 0, 0, 0);
  vec3 n = normalize(intersection->normal);
  point3 pos=intersection->position;
  material_s* mat=intersection->mat;
  Ray rayOmbre;
  rayOmbre.tmin=0;
  for(int i = 0 ; i<ln; i++){
    l=normalize(l_s[i]->position-pos);
    rayOmbre.orig=pos+acne_eps*l;
    rayOmbre.dir=l;
    rayOmbre.tmax=dot(l_s[i]->position-n,l_s[i]->position-pos);
    if(!intersectSceneOmbre(scene,&rayOmbre)) ret+=shade(n,v,l,l_s[i]->color,mat); // ajouter test fonction si res >1
  }
  return ret;
}*/
color3 calculShade(Scene* scene,Intersection* intersection, vec3 v , Lights l_s){
  int ln = l_s.size();
  point3 pos = intersection->position;
  Ray rayOmbre;
  rayOmbre.tmin=0;
  color3 res=color3(0.f);
  vec3 n = normalize(intersection->normal);
  material_s* m=intersection->mat;
  vec3 l;
  for(int i = 0 ; i < ln ; i++){
    l=normalize(l_s[i]->position-pos);
    rayOmbre.orig=pos + acne_eps*l;
    rayOmbre.dir=l;
    rayOmbre.tmax=distance(pos,l_s[i]->position);
    if(!intersectSceneOmbre(scene,&rayOmbre)){ 
      res += shade(n,v,l,l_s[i]->color,m);
      }
  }
  return res;
}

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, int depth) {
  color3 ret = color3(0, 0, 0);
  Intersection intersection;
  Ray newRay;
  int ricochet = 10;
  if(depth>9){
    return color3(0.f);
  }
  if(intersectScene(scene,ray,&intersection)){
    ret+=calculShade(scene,&intersection,normalize(-(ray->dir)),scene->lights);
  }
  else{
    return scene->skyColor;
  }
  if(ret[0]>=1 && ret[1]>= 1 && ret[2]>= 1){
    return color3(1.f);
  }
  if( depth<ricochet){
    vec3 reflectL = reflect(ray->dir,intersection.normal);
    rayInit(&newRay,intersection.position+acne_eps*reflectL,reflectL,0,10000);
    ret+= RDM_Fresnel(dot(newRay.dir,normalize(-(ray->dir)+newRay.dir)),1,intersection.mat->IOR)*trace_ray(scene,&newRay, tree, depth+1)* intersection.mat->specularColor;
    }
  else{
    ret=scene->skyColor;
  }

  return ret;
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing
  //! and kdtree initializaion
  float aspect = 1.f / scene->cam.aspect;

  KdTree *tree = NULL;


//! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f);   //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
                     aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x =
      (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;


  for (size_t j = 0; j < img->height; j++) {
    if (j != 0)
      printf("\033[A\r");
    float progress = (float)j / img->height * 100.f;
    printf("progress\t[");
    int cpt = 0;
    for (cpt = 0; cpt < progress; cpt += 5)
      printf(".");
    for (; cpt < 100; cpt += 5)
      printf(" ");
    printf("]\n");
#pragma omp parallel for
    for (size_t i = 0; i < img->width; i++) {
      color3 *ptr = getPixelPtr(img, i, j);
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                     float(i) * dx + float(j) * dy;

      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      *ptr = trace_ray(scene, &rx, tree, 0);

    }
  }
}
