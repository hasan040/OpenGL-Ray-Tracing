#include<windows.h>
#include<GL/glut.h>
#include<iostream>
#include<stack>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<vector>
#include<unistd.h>

#include "1405040_classes.h"
#include "bitmap_image.hpp"

using namespace std;

//double NEG_INF = -999999.0;
//double POS_INF = 999999.0;

double windowWidth = 500;
double windowHeight = 500;

double imageWidth;
double imageHeight;

int drawaxes;

double rdn_angle = (pi * 3)/180; // 3 degree,in rad format
double _distance = 3;
double ax_length = 800;

double viewAngle = (pi * 90)/180;   // 100 degree,in rad format



point camPos;

point u_vec;
point r_vec;
point l_vec;


int pxDimension;
int objCount;
int ltSourceCount;



//starting methods here



//vector <Shape *> objectList;

//vector <Light *> lightList;

void freeObjects()
{
    for(unsigned int i=0;i<objectList.size();i++){
        delete objectList[i];
    }

    for(unsigned int i=0;i<lightList.size();i++){
        delete lightList[i];
    }
}

void drawAxes()
{
	if(drawaxes==1)
	{
	    //glColor3f(1, 1, 1);
		glBegin(GL_LINES);{
		    glColor3f(1, 0, 0);
			glVertex3f( ax_length,0,0);
			glVertex3f(-ax_length,0,0);



			glColor3f(0, 0, 1);
			glVertex3f(0,0, ax_length);
			glVertex3f(0,0,-ax_length);

            glColor3f(0,1,0);
			glVertex3f(0,-ax_length,0);
			glVertex3f(0, ax_length,0);
		}
		glEnd();


	}
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
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
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);


		for(j=0;j<slices;j++)
		{
		    if(j & 1)
                glColor3f(0,0,0);
            else
                glColor3f(1,1,1);
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
}


void drawFloor()
{
     glColor3f(1,0,0);

     double x,y,a;

     x = -200;
     y = -200;
     a = 30;

     double pre_x = x;

     double flag = 1;

     for(int i=1;i<=30;i++){
        for(int j=1;j<=30;j++){
            glColor3f(flag,flag,flag);
            glBegin(GL_QUADS);
            {
                glVertex3f(x,y,0);
                glVertex3f(x+a,y,0);
                glVertex3f(x+a,y+a,0);
                glVertex3f(x,y+a,0);
            }
            glEnd();
            x = x + a;

            if(flag)
                flag = 0;
            else
                flag = 1;
        }

        x = pre_x;
        y = y+a;

        if(flag)
            flag = 0;
        else
            flag = 1;

     }







}




void drawSS()
{


    for(unsigned int i=0;i<objectList.size();i++){
        Shape * sp = objectList[i];
        sp->draw();
    }

    for(unsigned int i=0;i<lightList.size();i++){
        Light * lt = lightList[i];
        point _p1 = lt->light_pos;
        glColor3f(lt->color[0],lt->color[1],lt->color[2]);
        glPushMatrix();
        glRotatef(90,1,0,0);
        glBegin(GL_QUADS);
        glVertex3f(_p1.x,_p1.y,_p1.z);
        glVertex3f(_p1.x+10,_p1.y,_p1.z);
        glVertex3f(_p1.x+10,_p1.y+10,_p1.z);
        glVertex3f(_p1.x,_p1.y+10,_p1.z);
        glEnd();
        glPopMatrix();
    }
}

void capture()
{

    cout << "  CAPTURING.........."<<endl;


    bitmap_image image(pxDimension,pxDimension);




    double planeDistance = (windowHeight / 2) / tan(viewAngle / 2);

    point topLeft = camPos + l_vec * planeDistance - r_vec * (windowWidth / 2) + u_vec * (windowHeight / 2);

    double du = windowWidth / imageWidth;
    double dv = windowHeight / imageHeight;

    topLeft = topLeft + r_vec * (0.5*du) - u_vec * (0.5*dv);


    for(int i=0;i<imageWidth;i++){
        for(int j=0;j<imageHeight;j++){


            point curPixVec = topLeft + r_vec * (i - 0.0) * du - u_vec * (j - 0.0) * dv;


            point _dir = curPixVec-camPos;//actual direction for the rest of the shape except checkerboard

            _dir = _dir.normalize(); // unit vector for the ray direction

            Ray ray(camPos,_dir);



            double t_min = POS_INF;

            point dummy_color;

            Shape * targetShape = objectList[0];

            for(unsigned int k=0;k<objectList.size();k++){
                double t;

                Shape * dummyShape = objectList[k];
                point temp_color;
                if(dummyShape->didIntersect(ray,t,temp_color) != -1.0){
                    if(t_min > t && t > 0.0){
                        t_min = t;
                        dummy_color = temp_color;

                        targetShape = objectList[k];
                    }
                }

            }

            if(t_min != POS_INF){
                //image.set_pixel(i,j,255*expectedShape->color[0],255*expectedShape->color[1],255*expectedShape->color[2]);
                targetShape->illuminate(ray,dummy_color,t_min,0);
                image.set_pixel(i,j,255.0*dummy_color.x,255.0*dummy_color.y,255.0*dummy_color.z);
            }








            //image.set_pixel(j,i,255*expectedShape->color[0],255*expectedShape->color[1],255*expectedShape->color[2]);



            /*
            Shape * sp = objectList[0];
            double t;
            if(sp->didIntersect(ray,t) != -1){
                counter++;
                image.set_pixel(j,i,255*sp->color[0],255*sp->color[1],255*sp->color[2]);
            }

            */


        }
    }




    cout << "\n.........FINISHED COMPUTATION..............!!"<<endl;

    image.save_image("E:/GraphicsCode/raster/demo.bmp");

}


void op_rotate(point & fx,point & rt1,point & rt2,double _rdn)
{


    point cross1,cross2;

    cross1.x = fx.y * rt1.z - fx.z * rt1.y;
    cross1.y = fx.z * rt1.x - fx.x * rt1.z;
    cross1.z = fx.x * rt1.y - fx.y * rt1.x;

    cross2.x = fx.y * rt2.z - fx.z * rt2.y;
    cross2.y = fx.z * rt2.x - fx.x * rt2.z;
    cross2.z = fx.x * rt2.y - fx.y * rt2.x;

    rt1.x = rt1.x * cos(_rdn) + cross1.x * sin(_rdn);
    rt1.y = rt1.y * cos(_rdn) + cross1.y * sin(_rdn);
    rt1.z = rt1.z * cos(_rdn) + cross1.z * sin(_rdn);

    rt2.x = rt2.x * cos(_rdn) + cross2.x * sin(_rdn);
    rt2.y = rt2.y * cos(_rdn) + cross2.y * sin(_rdn);
    rt2.z = rt2.z * cos(_rdn) + cross2.z * sin(_rdn);

}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){


        //for rotation
        case '2':
            op_rotate(u_vec,l_vec,r_vec,-rdn_angle);
            break;

        case '1':
            op_rotate(u_vec,l_vec,r_vec,rdn_angle);
            break;

        case '3':
            op_rotate(r_vec,u_vec,l_vec,rdn_angle);
            break;
        case '4':
            op_rotate(r_vec,u_vec,l_vec,-rdn_angle);
            break;
        case '5':
            op_rotate(l_vec,u_vec,r_vec,rdn_angle);
            break;
        case '6':
            op_rotate(l_vec,u_vec,r_vec,-rdn_angle);
            break;
        case '0':
            capture();
            break;


		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			camPos.x -= _distance*l_vec.x;
			camPos.y -= _distance*l_vec.y;
			camPos.z -= _distance*l_vec.z;


			break;
		case GLUT_KEY_UP:		// up arrow key
			camPos.x += _distance*l_vec.x;
			camPos.y += _distance*l_vec.y;
			camPos.z += _distance*l_vec.z;

			break;

		case GLUT_KEY_RIGHT:
			camPos.x += _distance*r_vec.x;
			camPos.y += _distance*r_vec.y;
			camPos.z += _distance*r_vec.z;

			break;
		case GLUT_KEY_LEFT:
			camPos.x -= _distance*r_vec.x;
			camPos.y -= _distance*r_vec.y;
			camPos.z -= _distance*r_vec.z;

			break;

		case GLUT_KEY_PAGE_UP:
			camPos.x += _distance*u_vec.x;
			camPos.y += _distance*u_vec.y;
			camPos.z += _distance*u_vec.z;

			break;
		case GLUT_KEY_PAGE_DOWN:
			camPos.x -= _distance*u_vec.x;
			camPos.y -= _distance*u_vec.y;
			camPos.z -= _distance*u_vec.z;

			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	1,0,0);




	gluLookAt(camPos.x,camPos.y,camPos.z,
           camPos.x + l_vec.x,camPos.y + l_vec.y,camPos.z + l_vec.z,
           u_vec.x,u_vec.y,u_vec.z);






	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();



    drawSS();




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization

	drawaxes=1;
	//cameraHeight=150.0;
	//cameraAngle=1.0;
	//angle=0;

	//added components for camera position

	camPos.x = 150;
	camPos.y = 150;
	camPos.z = 80;//earlier 0

	u_vec.x = 0;
	u_vec.y = 0;
	u_vec.z = 1;

	r_vec.x = -1/sqrt(2);
	r_vec.y = 1/sqrt(2);
	r_vec.z = 0;

	l_vec.x = -1/sqrt(2);
	l_vec.y = -1/sqrt(2);
	l_vec.z = 0;






	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective((180 * viewAngle)/pi,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void loadData()
{
    cout << "Reading File."<<endl;

    ifstream file;

    file.open("E:/GraphicsCode/raster/scene.txt",ios::in);

    string strline;

    getline(file,strline);
    recurLevel = atoi(strline.c_str());

    getline(file,strline);
    pxDimension = atof(strline.c_str());
    imageWidth = pxDimension;
    imageHeight = pxDimension;

    getline(file,strline);



    getline(file,strline);
    objCount = atoi(strline.c_str());

    //cout << "Total object :"<<objCount<<endl;

    for(int i=0;i<objCount;i++){
        getline(file,strline);



        if(strline == "sphere"){


            getline(file,strline); //center
            point center;
            stringSeparator3(strline,center.x,center.y,center.z);


            getline(file,strline);//radius
            double radius;
            radius = atof(strline.c_str());

            Shape * shape = new Sphere(center,radius);

            getline(file,strline);//color
            double r_color,g_color,b_color;
            stringSeparator3(strline,r_color,g_color,b_color);
            shape->setColor(r_color,g_color,b_color);

            getline(file,strline);//co efficient
            double amb,diff,spec,rrc;
            stringSeparator4(strline,amb,diff,spec,rrc);
            shape->setCoEfficients(amb,diff,spec,rrc);


            getline(file,strline);//shine
            int shining = atoi(strline.c_str());
            shape->setShine(shining);

            objectList.push_back(shape);



            getline(file,strline);


        }
        else if(strline == "triangle"){

            getline(file,strline);//first point
            point a,b,c;
            stringSeparator3(strline,a.x,a.y,a.z);

            getline(file,strline);//second point
            stringSeparator3(strline,b.x,b.y,b.z);

            getline(file,strline);//third point
            stringSeparator3(strline,c.x,c.y,c.z);

            Shape *shape = new Triangle(a,b,c);//initialize triangle

            getline(file,strline);//set color
            double r_color,g_color,b_color;
            stringSeparator3(strline,r_color,g_color,b_color);
            shape->setColor(r_color,g_color,b_color);

            getline(file,strline);//set co efficients
            double amb,diff,spec,rrc;
            stringSeparator4(strline,amb,diff,spec,rrc);
            shape->setCoEfficients(amb,diff,spec,rrc);

            getline(file,strline);//set shine
            int shiningness;
            shiningness = atoi(strline.c_str());
            shape->setShine(shiningness);

            objectList.push_back(shape);

            getline(file,strline);

        }
        else if(strline == "general"){

            getline(file,strline);
            double v1,v2,v3,v4,v5,v6,v7,v8,v9,v10;
            stringSeparator10(strline,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10);




            getline(file,strline);//......
            double x1,x2,x3,x4,x5,x6;
            stringSeparator6(strline,x1,x2,x3,x4,x5,x6);
            point _ref(x1,x2,x3);
            point _dim(x4,x5,x6);

            Shape * shape = new General(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,_ref,_dim);


            getline(file,strline); //color input:3
            double r_color,g_color,b_color;
            stringSeparator3(strline,r_color,g_color,b_color);
            shape->setColor(r_color,g_color,b_color);

            getline(file,strline);//co efficients input 4
            double amb,diff,spec,rrc;
            stringSeparator4(strline,amb,diff,spec,rrc);
            shape->setCoEfficients(amb,diff,spec,rrc);

            getline(file,strline);//shine input
            int shiningness;
            shiningness = atoi(strline.c_str());
            shape->setShine(shiningness);


            objectList.push_back(shape);


            getline(file,strline);

        }
    }

    getline(file,strline);

    ltSourceCount = atoi(strline.c_str());

    //cout << "total light source :"<<ltSourceCount<<endl;

    for(int i=0;i<ltSourceCount;i++){
        getline(file,strline);
        double a,b,c;
        stringSeparator3(strline,a,b,c);

        point loc(a,b,c);


        getline(file,strline);
        double _r,_g,_b;
        stringSeparator3(strline,_r,_g,_b);

        Light * light = new Light(loc,_r,_g,_b);

        lightList.push_back(light);


    }

    //cout << "last light :"<<lightList[lightList.size()-1]->light_pos.z<<endl;

    point p1(-250.0,-250.0,0.0);

    //Shape * fr = new MyFloor(p1,30.0);

    Shape * fr = new Floor(2000,30.0);

    objectList.push_back(fr);

    cout << "Finish reading file.\n"<<endl;

    file.close();

}

int main(int argc, char **argv){

    loadData();

	glutInit(&argc,argv);

	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray Tracing");

	init();


	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	freeObjects();

	return 0;
}



