//#ifndef 1405040_CLASSES_H_INCLUDED
//#define 1405040_CLASSES_H_INCLUDED

#ifndef INCLUDE_1405040_CLASSES_H
#define INCLUDE_1405040_CLASSES_H
#include <iostream>


using namespace std;

#define pi (2*acos(0.0))

double NEG_INF = -999999.0;
double POS_INF = 999999.0;
int recurLevel;



struct point
{
	double x,y,z;

	point(){}

	point(double a,double b,double c){
	    x = a;
	    y = b;
	    z = c;
	}

    point operator+(const point & temp)const{
        point t;
        t.x = x + temp.x;
        t.y = y + temp.y;
        t.z = z + temp.z;

        return t;
    }

    point operator-(const point & temp)const{
        point t;
        t.x = x - temp.x;
        t.y = y - temp.y;
        t.z = z - temp.z;

        return t;
    }

    point normalize(){
        point t;
        double mgn = sqrt(x*x + y*y + z*z);
        t.x = x/mgn;
        t.y = y/mgn;
        t.z = z/mgn;
        return t;
    }

    double dot_product(const point & temp)const{
        return x * temp.x + y * temp.y + z * temp.z;
    }

    point operator*(const double & temp)const{
        point t;
        t.x = x * temp;
        t.y = y * temp;
        t.z = z * temp;
        return t;
    }

    point operator/(const double & temp)const{
        point t;
        t.x = x / temp;
        t.y = y / temp;
        t.z = z / temp;
        return t;
    }

    point cross_product(const point & temp)const{
        point t;
        t.x = y * temp.z - z * temp.y;
        t.y = z * temp.x - x * temp.z;
        t.z = x * temp.y - y * temp.x;
        return t;
    }

    point nonvector_product(const point & temp)const{
        point t;
        t.x = x * temp.x;
        t.y = y * temp.y;
        t.z = z * temp.z;
        return t;
    }


};


class MyQuad{
public:
    point a_point,b_point,c_point,d_point;
    MyQuad(point w,point x,point y,point z){
        a_point = w;
        b_point = x;
        c_point = y;
        d_point = z;
    }
};

class Light{
public:
    point light_pos;
    double color[3];
    Light(){}
    virtual ~Light(){}
    Light(point pos,double r,double g,double b){
        light_pos = pos;
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
};


class Ray{
public:
    point start;
    point direction;
    Ray(){}
    Ray(point init,point dir){
        start = init;
        direction = dir;//normalized form should be received
        //direction = direction.normalize();
    }

};

class Shape
{
public:
    point reference_point;
    point ref_point2;
    point ref_point3;
    double height,width,length;
    double color[3];
    double coEfficients[4];
    int shine;
    Shape(point p,double d){}
    Shape(point a,point b,point c){}
    Shape(double,double){}
    Shape(double _a1,double _b2,double _c3,double _d4,double _e5,double _f6,double _g7,double _h8,double _i9,double _j10,point _ref,point _dim){}
    virtual ~Shape(){}
    virtual void draw()=0;
    virtual point getNormalUnitVec(point) = 0;
    virtual double didIntersect(Ray & ,double &,point &) = 0;
    virtual void illuminate(Ray & ,point &, double & ,int ) = 0;
    void setColor(double r,double g,double b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
    void setShine(int shiningness){
        shine = shiningness;
    }
    void setCoEfficients(double amb,double diff,double spec,double rc){
        coEfficients[0] = amb;
        coEfficients[1] = diff;
        coEfficients[2] = spec;
        coEfficients[3] = rc;
    }


};


/**********************
/ GLOBAL VECTOR LIST
**********************/
vector <Shape *> objectList;

vector <Light *> lightList;

/**********************
/
**********************/

class Sphere : public Shape
{
public:

    Sphere(point center,double radius):Shape(center,radius){
        reference_point = center;
        length = radius;
        point demo_point(-1.0,-1.0,-1.0);
        ref_point2 = demo_point;
        ref_point3 = demo_point;

        height = radius;
        width = radius;
    }

    virtual void draw(){

        glPushMatrix();
        glTranslatef(reference_point.x,reference_point.y,reference_point.z);

	    struct point points[100][100];
	    int i,j;
	    int stacks,slices;
	    double h,r;

	    stacks = 70;
	    slices = 70;

	    //generate points
	    for(i=0;i<=stacks;i++)
	    {
		    h=length*sin(((double)i/(double)stacks)*(pi/2));
		    r=length*cos(((double)i/(double)stacks)*(pi/2));
		    for(j=0;j<=slices;j++)
		    {
			    points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			    points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			    points[i][j].z=h;
		    }
	     }
	//draw quads using generated points
	     for(i=0;i<stacks;i++)
	     {
             glColor3f(color[0],color[1],color[2]);


		     for(j=0;j<slices;j++)
             {
                 glBegin(GL_QUADS);{
                     //upper hemisphere
                     glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				     glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				     glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				     glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                     //lower hemisphere
                     glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				     glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				     glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				     glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			     }glEnd();
		      }
	     }

	     glPopMatrix();
    }

    virtual point getNormalUnitVec(point givenVec){
        point norm = givenVec - reference_point;
        return norm.normalize();

    }
    virtual double didIntersect(Ray & ray, double & t, point & dummy_color)
    {

        point R_dir = ray.direction;
        point R_org = ray.start - reference_point;

        double b = 2.0 * R_dir.dot_product(R_org);
        double c = R_org.dot_product(R_org) - length * length;
        double dis = b*b - 4.0*c;
        if (dis < 0.0)
            return -1.0;
        dis = sqrt(dis);
        double t1 = -b - dis;
        double t2 = -b + dis;

        double t_min,t_max;


        t_min = min(t1,t2);
        t_max = max(t1,t2);


        point _temp(color[0],color[1],color[2]);
        dummy_color = _temp;
        if(t_min >= 1e-4f){
            t = t_min/2.0;

            return t;
        }
        else if(t_max >= 1e-4f && t_min < 1e-4f){
            t = t_max/2.0;
            return t;
        }

        else{
            return -1.0;
        }
    }

    virtual void illuminate(Ray & ray,point & dummy_color,double & t_min, int level){
        dummy_color = dummy_color * coEfficients[0];

        if(level >= recurLevel){
            return;
        }

        point intersectionPoint = ray.start + ray.direction * t_min;

        point intersectionPointColor(color[0],color[1],color[2]);


        //computing for light sources

        point N_VEC = getNormalUnitVec(intersectionPoint);

        for(unsigned int i=0;i<lightList.size();i++){
            Light * light = lightList[i];

            //check if the light falls directly

            point light_direction = light->light_pos - intersectionPoint;
            light_direction = light_direction.normalize();
            Ray _light_ray(intersectionPoint,light_direction);

            double t_min_light = POS_INF;
            point _light_dummy_color;

            for(unsigned int k=0;k<objectList.size();k++){
                double _t;

                Shape * dummyShape = objectList[k];
                point _temp_color;
                if(dummyShape->didIntersect(_light_ray,_t,_temp_color) != -1.0){
                    if(t_min_light > _t && _t > 0.0){
                        t_min_light = _t;
                        _light_dummy_color = _temp_color;
                    }
                 }

             }

             if(t_min_light == POS_INF){

                point L_VEC = light->light_pos - intersectionPoint;
                L_VEC = L_VEC.normalize();
                double LambertValue = N_VEC.dot_product(L_VEC);

                point light_color(light->color[0],light->color[1],light->color[2]);

                dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * LambertValue * coEfficients[1];

                point R_VEC = N_VEC * (LambertValue * 2.0) - L_VEC;

                R_VEC = R_VEC.normalize();

                double phongValue = R_VEC.dot_product(ray.direction);

                double _phong = pow(phongValue,shine * 1.0);

                dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * _phong * coEfficients[2];
             }



            //end of the light tracing



        }

        //calculation for secondary ray

        point _distort(.001,.001,.001);

        point _mod = intersectionPoint + _distort;
        point _mod_dir = _mod - ray.start;
        _mod_dir = _mod_dir.normalize();

        double _temp_val = _mod_dir.dot_product(N_VEC);
        _temp_val = _temp_val * 2.0;

        point _reflected = _mod_dir - N_VEC * _temp_val;
        _reflected = _reflected.normalize();


        Ray ray_secondary(_mod,_reflected);




        double t_min_secondary = POS_INF;

        point dummy_color_secondary;



        for(unsigned int k=0;k<objectList.size();k++){
            double t;

            Shape * dummyShape = objectList[k];
            point temp_color;
            if(dummyShape->didIntersect(ray_secondary,t,temp_color) != -1.0){
                if(t_min_secondary > t && t > 0.0){
                    t_min_secondary = t;
                    dummy_color_secondary = temp_color;
                }
            }

        }

        if(t_min_secondary != POS_INF){

            dummy_color = dummy_color + dummy_color_secondary * coEfficients[3];
            illuminate(ray_secondary,dummy_color_secondary,t_min_secondary,level+1);


        }
        else{
            return;
        }

        //secondary ray pass
    }
};

class Triangle : public Shape
{
public:
    Triangle(point a,point b,point c):Shape(a,b,c){
        reference_point = a;
        ref_point2 = b;
        ref_point3 = c;
        height = -1;
        width = -1;
        length = -1;

    }

    virtual void draw(){

        glBegin(GL_TRIANGLES);
        {
            glColor3f(color[0],color[1],color[2]);
            glVertex3f(reference_point.x,reference_point.y,reference_point.z);
            glVertex3f(ref_point2.x,ref_point2.y,ref_point2.z);
            glVertex3f(ref_point3.x,ref_point3.y,ref_point3.z);

        }
        glEnd();

    }

    virtual point getNormalUnitVec(point givenVec)
    {
        point ba = ref_point2 - reference_point;
        point ca = ref_point3 - reference_point;

        point df = ba.cross_product(ca);
        df = df.normalize();

       // return (ba - ca).normalize();
        return df;

    }
    virtual double didIntersect(Ray & ray,double & t,point & dummy_color)
    {
        const double EPSILON = 1e-4f;

        point edge1,edge2,h,s,q;

        double a,f,u,v;

        edge1 = ref_point2 - reference_point;
        edge2 = ref_point3 - reference_point;

        h = ray.direction.cross_product(edge2);
        a = edge1.dot_product(h);

        if(a > -EPSILON && a < EPSILON)
            return -1.0; //parallel to the triangle
        f = 1.0 / a;
        s = ray.start - reference_point;
        u = f * s.dot_product(h);

        if(u < 0.0 || u > 1.0)
            return -1.0;
        q = s.cross_product(edge1);
        v = f * ray.direction.dot_product(q);

        if(v < 0.0 || (u + v) > 1.0)
            return -1.0;
        t = f * edge2.dot_product(q);

        if(t >= EPSILON){
            point _temp(color[0],color[1],color[2]);
            dummy_color = _temp;

            return t;
        }
        else{
            return -1.0;
        }

    }

    virtual void illuminate(Ray & ray,point & dummy_color,double & t_min, int  level){

        dummy_color = dummy_color * coEfficients[0];

        if(level >= recurLevel){
            return;
        }

        point intersectionPoint = ray.start + ray.direction * t_min;

        point intersectionPointColor(color[0],color[1],color[2]);



        point N_VEC = getNormalUnitVec(intersectionPoint);

        for(unsigned int i=0;i<lightList.size();i++){
            Light * light = lightList[i];

            //....start

            point light_direction = light->light_pos - intersectionPoint;
            light_direction = light_direction.normalize();
            Ray _light_ray(intersectionPoint,light_direction);

            double t_min_light = POS_INF;
            point _light_dummy_color;

            for(unsigned int k=0;k<objectList.size();k++){
                double _t;

                Shape * dummyShape = objectList[k];
                point _temp_color;
                if(dummyShape->didIntersect(_light_ray,_t,_temp_color) != -1.0){
                    if(t_min_light > _t && _t > 0.0){
                        t_min_light = _t;
                        _light_dummy_color = _temp_color;
                    }
                 }

             }

             if(t_min_light == POS_INF){

                point L_VEC = light->light_pos - intersectionPoint;
                L_VEC = L_VEC.normalize();
                double LambertValue = N_VEC.dot_product(L_VEC);

                point light_color(light->color[0],light->color[1],light->color[2]);

                dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * LambertValue * coEfficients[1];

                point R_VEC = N_VEC * (LambertValue * 2.0) - L_VEC;



                R_VEC = R_VEC.normalize();

                double phongValue = R_VEC.dot_product(ray.direction);

                double _phong = pow(phongValue,shine * 1.0);

                dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * _phong * coEfficients[2];
             }
            //end.........

            /*

            point L_VEC = light->light_pos - intersectionPoint;
            L_VEC = L_VEC.normalize();
            double LambertValue = N_VEC.dot_product(L_VEC);

            point light_color(light->color[0],light->color[1],light->color[2]);

            dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * LambertValue * coEfficients[1];

            point R_VEC = N_VEC * (LambertValue * 2) - L_VEC;
            R_VEC = R_VEC.normalize();

            double phongValue = R_VEC.dot_product(ray.direction);

            double _phong = pow(phongValue,shine * 1.0);

            dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * _phong * coEfficients[2];

            */
        }

        //calculation for secondary reflected ray


        point _distort(.001,.001,.001);

        point _mod = intersectionPoint + _distort;
        point _mod_dir = _mod - ray.start;
        _mod_dir = _mod_dir.normalize();

        double _temp_val = _mod_dir.dot_product(N_VEC);
        _temp_val = _temp_val * 2.0;

        point _reflected = _mod_dir - N_VEC * _temp_val;
        _reflected = _reflected.normalize();


        Ray ray_secondary(_mod,_reflected);




        double t_min_secondary = POS_INF;

        point dummy_color_secondary;



        for(unsigned int k=0;k<objectList.size();k++){
            double t;

            Shape * dummyShape = objectList[k];
            point temp_color;
            if(dummyShape->didIntersect(ray_secondary,t,temp_color) != -1.0){
                if(t_min_secondary > t && t > 0.0){
                    t_min_secondary = t;
                    dummy_color_secondary = temp_color;
                }
            }

        }

        if(t_min_secondary != POS_INF){
            dummy_color = dummy_color + dummy_color_secondary * coEfficients[3];
            illuminate(ray_secondary,dummy_color_secondary,t_min_secondary,level+1);

        }

        //finishing calculation

    }
};



class MyFloor :public Shape
{
private:

    vector <MyQuad> tileList;
    vector <point> colorList;


public:


    MyFloor(point startTile,double tileSize):Shape(startTile,tileSize){
        reference_point = startTile;
        height = tileSize;
        width = tileSize;
        length = tileSize;

        double x,y,a;
        x = reference_point.x;
        y = reference_point.y;

        a = tileSize;

        double flag = 1.0;

        for(int i=1;i<=16;i++){
            for(int j=1;j<=16;j++){
                point color(flag,flag,flag);

                point a_loc(x,y,0);
                point b_loc(x+a,y,0);
                point c_loc(x+a,y+a,0);
                point d_loc(x,y+a,0);

                MyQuad mq(a_loc,b_loc,c_loc,d_loc);
                tileList.push_back(mq);
                colorList.push_back(color);

                x = x + a;
                if(flag)
                    flag = 0.0;
                else
                    flag = 1.0;


            }

            x = reference_point.x;
            y = y+a;

            if(flag)
                flag = 0.0;
            else
                flag = 1.0;
        }

    }

    virtual void draw(){

        for(unsigned int i=0;i<colorList.size();i++){
            point cl = colorList[i];
            glColor3f(cl.x,cl.y,cl.z);
            MyQuad mq = tileList[i];

            glBegin(GL_QUADS);
            {
                glVertex3f(mq.a_point.x,mq.a_point.y,mq.a_point.z);
                glVertex3f(mq.b_point.x,mq.b_point.y,mq.b_point.z);
                glVertex3f(mq.c_point.x,mq.c_point.y,mq.c_point.z);
                glVertex3f(mq.d_point.x,mq.d_point.y,mq.d_point.z);

            }
            glEnd();

        }

    }

    virtual point getNormalUnitVec(point givenVec){
        point t(0,0,1.0);
        return t;
    }

    double quadInterSect(Ray & ray, MyQuad & mq, double &t){
        point S1 = mq.d_point;
        point S2 = mq.c_point;
        point S3 = mq.a_point;

        point ds_21 = S2 - S1;
        point ds_31 = S3 - S1;



        point normal(0.0,0.0,1.0);
        //point normal = ds_21.cross_product(ds_31);

        //normal = normal.normalize();

        double dot_norm = normal.dot_product(ray.direction);

        double dot_n1 = dot_norm;

        dot_n1 = dot_n1 < 0.0 ? -dot_n1 : dot_n1;

        if(dot_n1 < 1e-4f){ //dot_n1 < 1e-4f;changed
            return -1.0;
        }

        t = -normal.dot_product(ray.start - S1)/dot_norm;
        point M = ray.start + ray.direction * t;


        point d_MS1 = M - S1;
        double u = d_MS1.dot_product(ds_21);
        double v = d_MS1.dot_product(ds_31);

        if(u >= 0.0 && u <= ds_21.dot_product(ds_21) &&
           v >= 0.0 && v <= ds_31.dot_product(ds_31) ){
               return t;
           }
        else{

            return -1.0;
        }


    }

    virtual double didIntersect(Ray & ray, double &t,point & dummy_color){
        bool status = false;
        t = 999999.0;

        for(unsigned int i=0;i<tileList.size();i++){
            double t1;
            if(quadInterSect(ray,tileList[i],t1) != -1.0){
                if(t > t1){
                    t = t1;
                    status = true;

                    dummy_color = colorList[i];
                }
            }
        }

        if(status && t>=1e-4f){
            return t;
        }
        else{
            return -1.0;
        }



    }


    virtual void illuminate(Ray & ray,point & dummy_color,double & t_min, int  level){

        point intersectionPoint = ray.start + ray.direction * t_min;


        point _target_color = dummy_color;


        dummy_color = _target_color * coEfficients[0];

        if(level >= recurLevel){
            return;
        }



        point N_VEC = getNormalUnitVec(intersectionPoint);

        for(unsigned int i=0;i<lightList.size();i++){
            Light * light = lightList[i];

            //start

            point light_direction = light->light_pos - intersectionPoint;
            light_direction = light_direction.normalize();
            Ray _light_ray(intersectionPoint,light_direction);

            double t_min_light = POS_INF;
            point _light_dummy_color;

            for(unsigned int k=0;k<objectList.size();k++){
                double _t;

                Shape * dummyShape = objectList[k];
                point _temp_color;
                if(dummyShape->didIntersect(_light_ray,_t,_temp_color) != -1.0){
                    if(t_min_light > _t && _t > 0.0){
                        t_min_light = _t;
                        _light_dummy_color = _temp_color;
                    }
                 }

             }

             if(t_min_light == POS_INF){

                point L_VEC = light->light_pos - intersectionPoint;
                L_VEC = L_VEC.normalize();
                double LambertValue = N_VEC.dot_product(L_VEC);

                point light_color(light->color[0],light->color[1],light->color[2]);

                dummy_color = dummy_color + (light_color.nonvector_product(_target_color)) * LambertValue * coEfficients[1];

                point R_VEC = N_VEC * (LambertValue * 2) - L_VEC;
                R_VEC = R_VEC.normalize();

                double phongValue = R_VEC.dot_product(ray.direction);

                double _phong = pow(phongValue,shine * 1.0);

                dummy_color = dummy_color + (light_color.nonvector_product(_target_color)) * _phong * coEfficients[2];
             }

            //end........

            /*

            point L_VEC = light->light_pos - intersectionPoint;
            L_VEC = L_VEC.normalize();
            double LambertValue = N_VEC.dot_product(L_VEC);

            point light_color(light->color[0],light->color[1],light->color[2]);

            dummy_color = dummy_color + (light_color.nonvector_product(_target_color)) * LambertValue * coEfficients[1];

            point R_VEC = N_VEC * (LambertValue * 2) - L_VEC;
            R_VEC = R_VEC.normalize();

            double phongValue = R_VEC.dot_product(ray.direction);

            double _phong = pow(phongValue,shine * 1.0);

            dummy_color = dummy_color + (light_color.nonvector_product(_target_color)) * _phong * coEfficients[2];

            */

        }

        //calculation for secondary reflected ray


        point _distort(.001,.001,.001);

        point _mod = intersectionPoint + _distort;
        point _mod_dir = _mod - ray.start;
        _mod_dir = _mod_dir.normalize();

        double _temp_val = _mod_dir.dot_product(N_VEC);
        _temp_val = _temp_val * 2.0;

        point _reflected = _mod_dir - N_VEC * _temp_val;
        _reflected = _reflected.normalize();


        Ray ray_secondary(_mod,_reflected);




        double t_min_secondary = POS_INF;

        point dummy_color_secondary;



        for(unsigned int k=0;k<objectList.size();k++){
            double t;

            Shape * dummyShape = objectList[k];
            point temp_color;
            if(dummyShape->didIntersect(ray_secondary,t,temp_color) != -1.0){
                if(t_min_secondary > t && t > 0.0){
                    t_min_secondary = t;
                    dummy_color_secondary = temp_color;
                }
            }

        }

        if(t_min_secondary != POS_INF){
            dummy_color = dummy_color + dummy_color_secondary * coEfficients[3];
            illuminate(ray_secondary,dummy_color_secondary,t_min_secondary,level+1);

        }
        else{
            return;
        }

        //done calculation...



    }

};




class Floor : public Shape {
  public:
    int n;

    double floorWidth, tileWidth;

    Floor(double floorsize, double tilesize):Shape(floorsize,tilesize){

        floorWidth = 1000.0;
        tileWidth = 20.0;

        point df(-floorWidth / 4.0, -floorWidth / 4.0, 0.0);
        reference_point = df;
        n = floorWidth/tileWidth;

        coEfficients[0] = 0.4;
        coEfficients[1] = 0.2;
        coEfficients[2] = 0.1;
        coEfficients[3] = 0.3;

        shine = 5;


    }

    virtual void draw(){

        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                glBegin(GL_QUADS);


                {
                    glColor3f((i+j)&1,(i+j)&1,(i+j)&1);

                    glVertex3f(reference_point.x + tileWidth * i,reference_point.y + tileWidth * j, 0.0);
                    glVertex3f(reference_point.x + tileWidth * (i + 1),reference_point.y + tileWidth * j, 0.0);
                    glVertex3f(reference_point.x + tileWidth * (i + 1),reference_point.y + tileWidth * (j + 1), 0.0);
                    glVertex3f(reference_point.x + tileWidth * i,reference_point.y + tileWidth * (j + 1), 0.0);


                }

                glEnd();
            }
        }

    }

    virtual double didIntersect(Ray& ray,double & t,point & dummy_color){
        if (ray.direction.z == 0.0)
            return -1.0;
        t = -(ray.start.z / ray.direction.z);
        point _point = ray.start + ray.direction * t;
        int cell_i = (_point.x - reference_point.x) / tileWidth;
        int cell_j = (_point.y - reference_point.y) / tileWidth;

        if(cell_i < 0 || cell_i > n-1){
            return -1.0;
        }

        if(cell_j < 0 || cell_j > n-1){
            return -1.0;
        }

        point _demo((cell_i+cell_j)&1,(cell_i+cell_j)&1,(cell_i+cell_j)&1);
        dummy_color = _demo;
        if(t >= 1e-4f){
            return t;
        }
        else{
            return -1.0;
        }

    }

    virtual point getNormalUnitVec(point p){
        point t(0,0,1.0);
        return t; // always z-axis
    }

    virtual void illuminate(Ray & ray,point & dummy_color,double & t_min, int  level){

        point intersectionPoint = ray.start + ray.direction * t_min;
        int cell_i = (intersectionPoint.x - reference_point.x) / tileWidth;
        int cell_j = (intersectionPoint.y - reference_point.y) / tileWidth;

        if(cell_i < 0 || cell_i > n-1){
            return;
        }

        if(cell_j < 0 || cell_j > n-1){
            return ;
        }

        point _target_color((cell_i+cell_j)&1,(cell_i+cell_j)&1,(cell_i+cell_j)&1);


        dummy_color = dummy_color * coEfficients[0];//_target_color

        if(level >= recurLevel){
            return;
        }



        point N_VEC = getNormalUnitVec(intersectionPoint);

        for(unsigned int i=0;i<lightList.size();i++){
            Light * light = lightList[i];

            //start

            point light_direction = light->light_pos - intersectionPoint;
            light_direction = light_direction.normalize();
            Ray _light_ray(intersectionPoint,light_direction);

            double t_min_light = POS_INF;
            point _light_dummy_color;

            for(unsigned int k=0;k<objectList.size();k++){
                double _t;

                Shape * dummyShape = objectList[k];
                point _temp_color;
                if(dummyShape->didIntersect(_light_ray,_t,_temp_color) != -1.0){
                    if(t_min_light > _t && _t > 0.0){
                        t_min_light = _t;
                        _light_dummy_color = _temp_color;
                    }
                 }

             }

             if(t_min_light == POS_INF){

                point L_VEC = light->light_pos - intersectionPoint;
                L_VEC = L_VEC.normalize();
                double LambertValue = N_VEC.dot_product(L_VEC);

                point light_color(light->color[0],light->color[1],light->color[2]);

                dummy_color = dummy_color + (light_color.nonvector_product(_target_color)) * LambertValue * coEfficients[1];

                point R_VEC = N_VEC * (LambertValue * 2) - L_VEC;
                R_VEC = R_VEC.normalize();

                double phongValue = R_VEC.dot_product(ray.direction);

                double _phong = pow(phongValue,shine * 1.0);

                dummy_color = dummy_color + (light_color.nonvector_product(_target_color)) * _phong * coEfficients[2];
             }

            //end........

            /*

            point L_VEC = light->light_pos - intersectionPoint;
            L_VEC = L_VEC.normalize();
            double LambertValue = N_VEC.dot_product(L_VEC);

            point light_color(light->color[0],light->color[1],light->color[2]);

            dummy_color = dummy_color + (light_color.nonvector_product(_target_color)) * LambertValue * coEfficients[1];

            point R_VEC = N_VEC * (LambertValue * 2) - L_VEC;
            R_VEC = R_VEC.normalize();

            double phongValue = R_VEC.dot_product(ray.direction);

            double _phong = pow(phongValue,shine * 1.0);

            dummy_color = dummy_color + (light_color.nonvector_product(_target_color)) * _phong * coEfficients[2];

            */

        }

        //calculation for secondary reflected ray


        point _distort(.001,.001,.001);

        point _mod = intersectionPoint + _distort;
        point _mod_dir = _mod - ray.start;
        _mod_dir = _mod_dir.normalize();

        double _temp_val = _mod_dir.dot_product(N_VEC);
        _temp_val = _temp_val * 2.0;

        point _reflected = _mod_dir - N_VEC * _temp_val;
        _reflected = _reflected.normalize();


        Ray ray_secondary(_mod,_reflected);




        double t_min_secondary = POS_INF;

        point dummy_color_secondary;



        for(unsigned int k=0;k<objectList.size();k++){
            double t;

            Shape * dummyShape = objectList[k];
            point temp_color;
            if(dummyShape->didIntersect(ray_secondary,t,temp_color) != -1.0){
                if(t_min_secondary > t && t > 0.0){
                    t_min_secondary = t;
                    dummy_color_secondary = temp_color;
                }
            }

        }

        if(t_min_secondary != POS_INF){
            dummy_color = dummy_color + dummy_color_secondary * coEfficients[3];
            illuminate(ray_secondary,dummy_color_secondary,t_min_secondary,level+1);

        }
        else{
            return;
        }

        //done calculation...



    }

};



class General : public Shape {

private:
double A_1,B_2,C_3,D_4,E_5,F_6,G_7,H_8,I_9,J_10;

double lower_x,higher_x,lower_y,higher_y,lower_z,higher_z;


public:

    General(double _a1,double _b2,double _c3,double _d4,double _e5,double _f6,double _g7,double _h8,double _i9,double _j10,point _ref,point _dim):
        Shape(_a1,_b2,_c3,_d4,_e5,_f6,_g7,_h8,_i9,_j10,_ref,_dim){
            A_1 = _a1;
            B_2 = _b2;
            C_3 = _c3;
            D_4 = _d4;
            E_5 = _e5;
            F_6 = _f6;
            G_7 = _g7;
            H_8 = _h8;
            I_9 = _i9;
            J_10 = _j10;

            reference_point = _ref;
            length = _dim.x;
            width = _dim.y;
            height = _dim.z;

            //......

            lower_x = min(reference_point.x,reference_point.x+length);
            higher_x = max(reference_point.x,reference_point.x+length);

            lower_y = min(reference_point.y,reference_point.y+width);
            higher_y = max(reference_point.y,reference_point.y+width);

            lower_z = min(reference_point.z,reference_point.z+height);
            higher_z = max(reference_point.z,reference_point.z+height);

        }

    virtual void draw(){

    }

    bool inTheBox(point &p){

        int counter = 0;
        if(length == 0.0){
            counter++;
        }
        else{
            if(p.x>=lower_x && p.x<=higher_x){
                counter++;
            }
        }
        if(width == 0.0){
            counter++;
        }
        else{
            if(p.y>=lower_y && p.y<=higher_y){
                counter++;
            }
        }
        if(height == 0.0){
            counter++;
        }
        else{
            if(p.z>=lower_z && p.z<=higher_z){
                counter++;
            }
        }

        if(counter == 3){
            return true;
        }
        return false;

    }





    virtual double didIntersect(Ray& ray,double & t,point & dummy_color){

        double xo = ray.start.x;
        double yo = ray.start.y;
        double zo = ray.start.z;

        double xd = ray.direction.x;
        double yd = ray.direction.y;
        double zd = ray.direction.z;

        double aq = A_1*xd*xd + B_2*yd*yd + C_3*zd*zd + D_4*xd*yd + E_5*xd*zd + F_6*yd*zd;
        double bq = 2.0*A_1*xo*xd + 2.0*B_2*yo*yd + 2.0*C_3*zo*zd + D_4*(xo*yd + yo*xd) + E_5*(xo*zd + zo*xd) + F_6*(yo*zd + yd*zo) + G_7*xd + H_8*yd + I_9*zd;
        double cq = A_1*xo*xo + B_2*yo*yo + C_3*zo*zo + D_4*xo*yo + E_5*xo*zo + F_6*yo*zo + G_7*xo + H_8*yo + I_9*zo + J_10;

        double dis = bq*bq - 4.0*aq*cq;

        if(dis < 0.0)
            return -1.0;

        dis = sqrt(dis);

        double t1 = -bq + dis;
        double t2 = -bq - dis;

        double _t_min = min(t1,t2)/(2.0*aq);
        double _t_max = max(t1,t2)/(2.0*aq);

        point _temp(color[0],color[1],color[2]);
        dummy_color = _temp;


        point INS_POINT;



        if(_t_min>=1e-4 ){
            INS_POINT = ray.start + ray.direction * _t_min;
            if(inTheBox(INS_POINT)){
                t = _t_min;
                return t;
            }
        }
        if(_t_max >= 1e-4 ){
            INS_POINT = ray.start + ray.direction * _t_max;
            if(inTheBox(INS_POINT)){
                t = _t_max;
                return t;
            }

        }

        return -1.0;
    }

    virtual point getNormalUnitVec(point p){
        point t;
        t.x = 2*A_1*p.x + D_4*p.y + E_5*p.z + G_7;
        t.y = 2*B_2*p.y + D_4*p.x + F_6*p.z + H_8;
        t.z = 2*C_3*p.z + E_5*p.x + F_6*p.y + I_9;
        t = t.normalize();
        return t; // always z-axis
    }

    virtual void illuminate(Ray & ray,point & dummy_color,double & t_min, int  level){
        dummy_color = dummy_color * coEfficients[0];

        if(level >= recurLevel){
            return;
        }

        point intersectionPoint = ray.start + ray.direction * t_min;

        point intersectionPointColor(color[0],color[1],color[2]);


        point N_VEC = getNormalUnitVec(intersectionPoint);

        for(unsigned int i=0;i<lightList.size();i++){
            Light * light = lightList[i];

            //start.....

            point light_direction = light->light_pos - intersectionPoint;
            light_direction = light_direction.normalize();
            Ray _light_ray(intersectionPoint,light_direction);

            double t_min_light = POS_INF;
            point _light_dummy_color;

            for(unsigned int k=0;k<objectList.size();k++){
                double _t;

                Shape * dummyShape = objectList[k];
                point _temp_color;
                if(dummyShape->didIntersect(_light_ray,_t,_temp_color) != -1.0){
                    if(t_min_light > _t && _t > 0.0){
                        t_min_light = _t;
                        _light_dummy_color = _temp_color;
                    }
                 }

             }

             if(t_min_light == POS_INF){

                point L_VEC = light->light_pos - intersectionPoint;
                L_VEC = L_VEC.normalize();
                double LambertValue = N_VEC.dot_product(L_VEC);

                point light_color(light->color[0],light->color[1],light->color[2]);

                dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * LambertValue * coEfficients[1];

                point R_VEC = N_VEC * (LambertValue * 2.0) - L_VEC;

                R_VEC = R_VEC.normalize();

                double phongValue = R_VEC.dot_product(ray.direction);

                double _phong = pow(phongValue,shine * 1.0);

                dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * _phong * coEfficients[2];
             }


            //end.....



            /*

            point L_VEC = light->light_pos - intersectionPoint;
            L_VEC = L_VEC.normalize();
            double LambertValue = N_VEC.dot_product(L_VEC);

            point light_color(light->color[0],light->color[1],light->color[2]);

            dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * LambertValue * coEfficients[1];

            point R_VEC = N_VEC * (LambertValue * 2) - L_VEC;
            R_VEC = R_VEC.normalize();

            double phongValue = R_VEC.dot_product(ray.direction);

            double _phong = pow(phongValue,shine * 1.0);

            dummy_color = dummy_color + (light_color.nonvector_product(intersectionPointColor)) * _phong * coEfficients[2];

            */

        }

        //secondary reflected ray calculation starts here


        point _distort(.001,.001,.001);

        point _mod = intersectionPoint + _distort;
        point _mod_dir = _mod - ray.start;
        _mod_dir = _mod_dir.normalize();

        double _temp_val = _mod_dir.dot_product(N_VEC);
        _temp_val = _temp_val * 2.0;

        point _reflected = _mod_dir - N_VEC * _temp_val;
        _reflected = _reflected.normalize();


        Ray ray_secondary(_mod,_reflected);




        double t_min_secondary = POS_INF;

        point dummy_color_secondary;



        for(unsigned int k=0;k<objectList.size();k++){
            double t;

            Shape * dummyShape = objectList[k];
            point temp_color;
            if(dummyShape->didIntersect(ray_secondary,t,temp_color) != -1.0){
                if(t_min_secondary > t && t > 0.0){
                    t_min_secondary = t;
                    dummy_color_secondary = temp_color;
                }
            }

        }

        if(t_min_secondary != POS_INF){
            dummy_color = dummy_color + dummy_color_secondary * coEfficients[3];
            illuminate(ray_secondary,dummy_color_secondary,t_min_secondary,level+1);

        }
        else{
            return;
        }

        //end calculation of the secondary ray

    }


};

//

//



void stringSeparator2(string strline,double & v1,double & v2)
{
    string arr[2];

    int i;
    int strlength = strline.size();
    int index = 0;

    string temp = "";

    for(i=0;i<strlength;i++){
        if(strline[i] == ' '){
            arr[index] = temp;
            index++;
            temp = "";

        }
        else{
            temp += strline[i];
        }
    }
    arr[index] = temp;

    v1 = atof(arr[0].c_str());
    v2 = atof(arr[1].c_str());

}

void stringSeparator3(string strline,double & v1,double & v2,double & v3)
{
    string arr[3];

    int i;
    int strlength = strline.size();
    int index = 0;

    string temp = "";

    for(i=0;i<strlength;i++){
        if(strline[i] == ' '){
            arr[index] = temp;
            index++;
            temp = "";

        }
        else{
            temp += strline[i];
        }
    }
    arr[index] = temp;

    v1 = atof(arr[0].c_str());
    v2 = atof(arr[1].c_str());
    v3 = atof(arr[2].c_str());


}

void stringSeparator4(string strline,double & v1,double & v2,double & v3,double & v4)
{

    string arr[4];

    int i;
    int strlength = strline.size();
    int index = 0;

    string temp = "";

    for(i=0;i<strlength;i++){
        if(strline[i] == ' '){
            arr[index] = temp;
            index++;
            temp = "";

        }
        else{
            temp += strline[i];
        }
    }
    arr[index] = temp;

    v1 = atof(arr[0].c_str());
    v2 = atof(arr[1].c_str());
    v3 = atof(arr[2].c_str());
    v4 = atof(arr[3].c_str());

}


void stringSeparator6(string strline,double & v1,double & v2,double & v3,double & v4,double & v5,double & v6)
{
    string arr[6];

    int i;
    int strlength = strline.size();
    int index = 0;

    string temp = "";

    for(i=0;i<strlength;i++){
        if(strline[i] == ' '){
            arr[index] = temp;
            index++;
            temp = "";

        }
        else{
            temp += strline[i];
        }
    }
    arr[index] = temp;

    v1 = atof(arr[0].c_str());
    v2 = atof(arr[1].c_str());
    v3 = atof(arr[2].c_str());
    v4 = atof(arr[3].c_str());

    v5 = atof(arr[4].c_str());
    v6 = atof(arr[5].c_str());

}


void stringSeparator10(string strline,double & v1,double & v2,double & v3,double & v4,double & v5,double & v6,double & v7,double & v8,double & v9,double & v10)
{
    string arr[10];

    int i;
    int strlength = strline.size();
    int index = 0;

    string temp = "";

    for(i=0;i<strlength;i++){
        if(strline[i] == ' '){
            arr[index] = temp;
            index++;
            temp = "";

        }
        else{
            temp += strline[i];
        }
    }
    arr[index] = temp;

    v1 = atof(arr[0].c_str());
    v2 = atof(arr[1].c_str());
    v3 = atof(arr[2].c_str());
    v4 = atof(arr[3].c_str());

    v5 = atof(arr[4].c_str());
    v6 = atof(arr[5].c_str());
    v7 = atof(arr[6].c_str());
    v8 = atof(arr[7].c_str());

    v9 = atof(arr[8].c_str());
    v10 = atof(arr[9].c_str());


}


#endif
