#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <mmsystem.h>
#include <math.h>
#include <cmath>
#include <bits/stdc++.h>
#include "include/BmpLoader.h"
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <sstream>
#include <vector>
void displayClock();
using namespace std;

int hours = 120;
int minutes = 0;
int seconds = 0;
float laps =0;


#define PI 3.14159265389
//for brazier curve
int anglex= 0, angley = 0, anglez = 0;
int wired=0;
int shcpt=1;
int animat = 0;
const int L=20;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;	//number of slices along x-direction
const int ntheta = 20;
const int width = 800;
const int height = 600;
//const float rat = 1.4*width/height;
const float rat = 1.4*width/height;

GLfloat eyeX = 2.5;
GLfloat eyeY = 2.5;
GLfloat eyeZ = 20;

GLfloat lookX = 1;
GLfloat lookY = -0.5;
GLfloat lookZ = -11;

float rot = 0;
int wndw_rot=0;
float inter_wndw_rot=0;
//float boat_z = 7.5;
float boat_z = 14.5;
float boat_rotate = 0;
float boat_x = 2;
float boat_y = -2.5;

bool light1=true, light2=false, light3=false;
bool diff_light_on=true, spec_light_on=false, amb_light_on=true;
bool target_on=true,fan_on=true;
unsigned int sceneCall = 0;
unsigned int ID[] = {};
unsigned int brick, grass;
float rotation_test = 0;
bool cat_rotate= false;
double hx=0,hy=1,hz=0;
float wcsClkDn[3],wcsClkUp[3];

static void resize(int width, int height)
{

const float ar = (float) width / (float) height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-8.0, 8.0, -8.0*(GLfloat)height/(GLfloat)width, 8.0*(GLfloat)height/(GLfloat)width, 2.0, 25.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;
    gluLookAt(eyeX,eyeY,eyeZ,lookX,lookY,lookZ,hx,hy,hz);

}

GLfloat ctrlpoints[L+1][3] =
{
    { 0.0, 0.0, 0.0}, { -0.3, 0.5, 0.0},
    { 0.1, 1.7, 0.0},{ 0.5, 1.5, 0.0},
    {1.0, 1.5, 0.0}, {1.4, 1.4, 0.0},
    {1.8, 0.4, 0.0},{2.2, 0.4, 0.0},
    {2.6, 1.5, 0.0}, {3.0, 1.4, 0.0},
    {3.4, 1.4, 0.0},{3.8, 1.4, 0.0},
    {4.2, 1.0, 0.0},{4.6, 1.0, 0.0},
    {5.0, 1.0, 0.0},{5.4, 1.0, 0.0},
    {5.8, 0.5, 0.0},{6.2, 0.5, 0.0},
    {6.6, 0.5, 0.0},{7.2, 0.2, 0.0},
    {6.8, 0.52, 0.0}
};




class point1
{
public:
    point1()
    {
        x=0;
        y=0;
    }
    int x;
    int y;
} clkpt[2];



void scsToWcs(float sx,float sy, float wcsv[3] );
void processMouse(int button, int state, int x, int y);
int flag=0;
GLint viewport[4]; //var to hold the viewport info
GLdouble modelview[16]; //var to hold the modelview info
GLdouble projection[16]; //var to hold the projection matrix info
void scsToWcs(float sx,float sy, float wcsv[3] )
{

    GLfloat winX, winY, winZ; //variables to hold screen x,y,z coordinates
    GLdouble worldX, worldY, worldZ; //variables to hold world x,y,z coordinates

    glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
    glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
    glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

    winX = sx;
    winY = (float)viewport[3] - (float)sy;
    winZ = 0;

    //get the world coordinates from the screen coordinates
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);
    wcsv[0]=worldX;
    wcsv[1]=worldY;
    wcsv[2]=worldZ;


}
void processMouse(int button, int state, int x, int y)
{
    if(button==GLUT_LEFT_BUTTON && state==GLUT_DOWN)
    {
        if(flag!=1)
        {
            flag=1;
            clkpt[0].x=x;
            clkpt[0].y=y;
        }


        scsToWcs(clkpt[0].x,clkpt[0].y,wcsClkDn);
        //cout<<"\nD: "<<x<<" "<<y<<" wcs: "<<wcsClkDn[0]<<" "<<wcsClkDn[1];
        cout<<"\nD: "<<x<<" "<<y<<" wcs: "<<wcsClkDn[0]<<" "<<wcsClkDn[1];
    }
    else if(button==GLUT_LEFT_BUTTON && state==GLUT_UP)
    {
        if (flag==1)
        {
            clkpt[1].x=x;
            clkpt[1].y=y;
            flag=0;
        }
        float wcs[3];
        scsToWcs(clkpt[1].x,clkpt[1].y,wcsClkUp);
        cout<<"\nU: "<<x<<" "<<y<<" wcs: "<<wcsClkUp[0]<<" "<<wcsClkUp[1];

        clikd=!clikd;
    }
}

//control points
long long nCr(int n, int r)
{
    if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}
void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}
void bottleBezier(int ntnew)
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < ntnew; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa ;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}

void insectBezier(int ntnew)
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < 10; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa + r*sina ;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);


        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}






static GLfloat v_cube[8][3] =
{
    {0,0,0},
    {0,0,1},
    {0,1,0},
    {0,1,1},

    {1,0,0},
    {1,0,1},
    {1,1,0},
    {1,1,1}
};

static GLubyte c_ind[6][4] =
{
    {3,1,5,7},  //front
    {6,4,0,2},  //back
    {2,3,7,6},  //top
    {1,0,4,5},  //bottom
    {7,5,4,6},  //right
    {2,0,1,3}   //left
};

static void getNormal3p(GLfloat x1, GLfloat y1, GLfloat z1,
                        GLfloat x2, GLfloat y2, GLfloat z2,
                        GLfloat x3, GLfloat y3, GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

void set_mat_prop(float colR, float colG, float colB, bool em=false, float shine=128)
{
    GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { colR, colG, colB, 1.0 };
    GLfloat mat_diffuse[] = { colR, colG, colB, 1.0 };
    GLfloat mat_specular[] = { colR, colG, colB, 1.0 };
    GLfloat mat_emission[] = {colR, colG, colB, 1.0};
    GLfloat mat_shininess[] = {shine};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    if(em)
        glMaterialfv( GL_FRONT, GL_EMISSION, mat_emission);
    else
        glMaterialfv( GL_FRONT, GL_EMISSION, no_mat);
}

void cube(float colR=0.5, float colG=0.5, float colB=0.5,
          bool em=false, float shine=128)
{
    set_mat_prop(colR,colG,colB,em,shine);

    glBegin(GL_QUADS);
    for (GLint i = 0; i <6; i++)
    {
        getNormal3p(v_cube[c_ind[i][0]][0], v_cube[c_ind[i][0]][1], v_cube[c_ind[i][0]][2],
                    v_cube[c_ind[i][1]][0], v_cube[c_ind[i][1]][1], v_cube[c_ind[i][1]][2],
                    v_cube[c_ind[i][2]][0], v_cube[c_ind[i][2]][1], v_cube[c_ind[i][2]][2]);

        glTexCoord2f(1,1);
        glVertex3fv(&v_cube[c_ind[i][0]][0]);
        glTexCoord2f(1,0);
        glVertex3fv(&v_cube[c_ind[i][1]][0]);
        glTexCoord2f(0,0);
        glVertex3fv(&v_cube[c_ind[i][2]][0]);
        glTexCoord2f(0,1);
        glVertex3fv(&v_cube[c_ind[i][3]][0]);
    }
    glEnd();
}

//pyramid
static GLfloat v_pyramid[5][3] =
{
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 2.0},
    {2.0, 0.0, 2.0},
    {2.0, 0.0, 0.0},
    {1.0, 4.0, 1.0}
};

static GLubyte p_Indices[4][3] =
{
    {4, 1, 2},
    {4, 2, 3},
    {4, 3, 0},
    {4, 0, 1}
};

static GLubyte quadIndices[1][4] =
{
    {0, 3, 2, 1}
};



void drawpyramid()
{

    glBegin(GL_TRIANGLES);
    for (GLint i = 0; i <4; i++)
    {
        //glColor3f(colors[i][0],colors[i][1],colors[i][2]);
        getNormal3p(v_pyramid[p_Indices[i][0]][0], v_pyramid[p_Indices[i][0]][1], v_pyramid[p_Indices[i][0]][2],
                    v_pyramid[p_Indices[i][1]][0], v_pyramid[p_Indices[i][1]][1], v_pyramid[p_Indices[i][1]][2],
                    v_pyramid[p_Indices[i][2]][0], v_pyramid[p_Indices[i][2]][1], v_pyramid[p_Indices[i][2]][2]);

        glVertex3fv(&v_pyramid[p_Indices[i][0]][0]);
        glVertex3fv(&v_pyramid[p_Indices[i][1]][0]);
        glVertex3fv(&v_pyramid[p_Indices[i][2]][0]);
    }
    glEnd();

    glBegin(GL_QUADS);
    for (GLint i = 0; i <1; i++)
    {
        //glColor3f(colors[4][0],colors[4][1],colors[4][2]);
        getNormal3p(v_pyramid[quadIndices[i][0]][0], v_pyramid[quadIndices[i][0]][1], v_pyramid[quadIndices[i][0]][2],
                    v_pyramid[quadIndices[i][1]][0], v_pyramid[quadIndices[i][1]][1], v_pyramid[quadIndices[i][1]][2],
                    v_pyramid[quadIndices[i][2]][0], v_pyramid[quadIndices[i][2]][1], v_pyramid[quadIndices[i][2]][2]);

        glVertex3fv(&v_pyramid[quadIndices[i][0]][0]);
        glVertex3fv(&v_pyramid[quadIndices[i][1]][0]);
        glVertex3fv(&v_pyramid[quadIndices[i][2]][0]);
        glVertex3fv(&v_pyramid[quadIndices[i][3]][0]);
    }
    glEnd();
    //glutSolidSphere (3.0, 20, 16);

}


//#define PI 3.1415927
void draw_cylinder(GLfloat radius,
                   GLfloat height,float R=0.0, float G=0, float B=0 )
{
    GLfloat x              = 0.0;
    GLfloat y              = 0.0;
    GLfloat angle          = 0.0;
    GLfloat angle_stepsize = 0.1;
    set_mat_prop(R,G,B);

    //tube
    //glColor3ub(R-40,G-40,B-40);
    glBegin(GL_QUAD_STRIP);
    angle = 0.0;
        while( angle < 2*PI ) {
            x = radius * cos(angle);
            y = radius * sin(angle);
            glVertex3f(x, y , height);
            glVertex3f(x, y , 0.0);
            angle = angle + angle_stepsize;
        }
        glVertex3f(radius, 0.0, height);
        glVertex3f(radius, 0.0, 0.0);
    glEnd();

    //circle
    //glColor3ub(R,G,B);
    glBegin(GL_POLYGON);
    angle = 0.0;
        while( angle < 2*PI ) {
            x = radius * cos(angle);
            y = radius * sin(angle);
            glVertex3f(x, y , height);
            angle = angle + angle_stepsize;
        }
        glVertex3f(radius, 0.0, height);
    glEnd();
}



void sphere(float rad, int slices, int stacks,
            float colR=0.5, float colG=0.5, float colB=0.5, bool em=false, float shine=128)
{
    set_mat_prop(colR,colG,colB,em,shine);
    glutSolidSphere(rad, slices, stacks);
}

bool right_light = false , activate_amb = true, activate_spec = true, activate_diff = true;

void light_right(float x_pos=1, float y_pos=-1, float z_pos=0)
{
    GLfloat no_light[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat light_amb[]  = {0.3, 0.3,0.3, 1.0};
    GLfloat light_diff[]  = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_spec[] = { 1.0, 0.0, 0.0, 1.0 };
    GLfloat light_pos[] = { x_pos, y_pos, z_pos, 1.0 };
    GLfloat light_attenuation[] = {1,0,0.002};


    if(right_light)
       glEnable( GL_LIGHT1);
    else
        glDisable(GL_LIGHT1);


        if(activate_amb)
            glLightfv( GL_LIGHT1, GL_AMBIENT, light_amb);
        else
            glLightfv( GL_LIGHT1, GL_AMBIENT, no_light);


        if(activate_diff)
            glLightfv( GL_LIGHT1, GL_DIFFUSE, light_diff);
        else
             glLightfv( GL_LIGHT1, GL_DIFFUSE, no_light);

        if(activate_spec)
            glLightfv( GL_LIGHT1, GL_SPECULAR, light_spec);
        else
            glLightfv( GL_LIGHT1, GL_SPECULAR, no_light);


    glLightfv( GL_LIGHT1, GL_POSITION, light_pos);
}
bool left_light = false;
void light_left(float x_pos, float y_pos, float z_pos)
{
    GLfloat no_light[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat light_amb[]  = {0.3, 0.3, 0.3, 1.0};
    GLfloat light_diff[]  = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_spec[] = { 1.0, 0.0, 0.0, 1.0 };
    GLfloat light_pos[] = { x_pos, y_pos, z_pos, 1.0 };
    GLfloat light_attenuation[] = {1,0,0.002};

       if(left_light)
         glEnable( GL_LIGHT2);
       else
          glDisable(GL_LIGHT2);

        if(activate_amb)
            glLightfv( GL_LIGHT2, GL_AMBIENT, light_amb);
        else
            glLightfv( GL_LIGHT2, GL_AMBIENT, no_light);


        if(activate_diff)
            glLightfv( GL_LIGHT2, GL_DIFFUSE, light_diff);
        else
             glLightfv( GL_LIGHT2, GL_DIFFUSE, no_light);

        if(activate_spec)
            glLightfv( GL_LIGHT2, GL_SPECULAR, light_spec);
        else
            glLightfv( GL_LIGHT2, GL_SPECULAR, no_light);


    glLightfv( GL_LIGHT2, GL_POSITION, light_pos);
}

bool spot_light = false;
void light_spot(float x_pos, float y_pos, float z_pos)
{
    GLfloat no_light[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat light_amb[]  = {0.3, 0.0, 0.0, 1.0};
    GLfloat light_diff[]  = { 1.0, 0.0, 0.0, 1.0 };
    GLfloat light_spec[] = { 1.0, 0.0, 0.0, 1.0 };
    GLfloat light_pos[] = { boat_x, -4, boat_z, 1.0 };
    GLfloat light_attenuation[] = {1,0,0.002};

       if(spot_light)
         glEnable( GL_LIGHT3);
       else
          glDisable(GL_LIGHT3);

        if(activate_amb)
            glLightfv( GL_LIGHT3, GL_AMBIENT, light_amb);
        else
            glLightfv( GL_LIGHT3, GL_AMBIENT, no_light);


        if(activate_diff)
            glLightfv( GL_LIGHT3, GL_DIFFUSE, light_diff);
        else
             glLightfv( GL_LIGHT3, GL_DIFFUSE, no_light);

        if(activate_spec)
            glLightfv( GL_LIGHT3, GL_SPECULAR, light_spec);
        else
            glLightfv( GL_LIGHT3, GL_SPECULAR, no_light);


    glLightfv( GL_LIGHT3, GL_POSITION, light_pos);

    GLfloat spot_direction[] = { boat_x-3, -5, 0 };
    glLightfv(GL_LIGHT3, GL_SPOT_DIRECTION, spot_direction);
    glLightf( GL_LIGHT3, GL_SPOT_CUTOFF, 90.0);
}


GLubyte* make_bw_tiles_texture(int tile_width, int tile_height, int tex_width=2048, int tex_height=2048)
{
    GLubyte* bw_tiles = new GLubyte[tex_width*tex_height*3];

    for(int i=0; i<tex_height; i++)
    {
        for(int j=0; j<tex_width; j++)
        {
            int c=(i/tile_width + j/tile_width) % 2;

            if(c==0)
            {
                bw_tiles[(i*tex_width+j)*3+0]=255;
                bw_tiles[(i*tex_width+j)*3+1]=255;
                bw_tiles[(i*tex_width+j)*3+2]=255;
            }
            else
            {
                bw_tiles[(i*tex_width+j)*3+0]=0;
                bw_tiles[(i*tex_width+j)*3+1]=0;
                bw_tiles[(i*tex_width+j)*3+2]=0;
            }
        }
    }

    return bw_tiles;
}

void LoadTexture(const char*filename,int num)
{
    glGenTextures(1, &ID[num]);
    glBindTexture(GL_TEXTURE_2D,  ID[num]);
    glPixelStorei(GL_UNPACK_ALIGNMENT,  ID[num]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}

void axes()
{
    float length = 15;
    float width = 0.3;

    // X-axis
    glPushMatrix();
    glTranslatef(length/2,0,0);
    glScalef(length,width,width);
    glTranslatef(-0.5,-0.5,-0.5);
    cube(0.8,0.1,0.1);
    glPopMatrix();

    // Y-axis
    glPushMatrix();
    glTranslatef(0,length/2,0);
    glScalef(width,length,width);
    glTranslatef(-0.5,-0.5,-0.5);
    cube(0.1,0.8,0.1);
    glPopMatrix();

    // Z-axis
    glPushMatrix();
    glTranslatef(0,0,length/2);
    glScalef(width,width,length);
    glTranslatef(-0.5,-0.5,-0.5);
    cube(0.1,0.1,0.8);
    glPopMatrix();
}


void showControlPoints()
{
    glPointSize(5.0);
    glColor3f(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i <=L; i++)
        glVertex3fv(&ctrlpoints[i][0]);
    glEnd();
}



void drawString(float x, float y, float z, string str) {
  glRasterPos3f(x, y, z);

  const char *cstr = str.c_str();

  for (const char* c = cstr; *c != '\0'; c++) {
        //char currentFont = *c;
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *c);
  }
}

int GAP =1000,j=0;
void timer(int seconds)
{
    if(j==0)
 {

    displayClock();
    minutes = minutes + 1;
    if(minutes == 1000)
        {hours-=1;
        minutes =0;

        }
    if(hours == 0)
        sceneCall=2;


 }

 else if(j==1)
 {
  minutes = minutes + 0;
 }
 else
 {
  minutes = 0;
 }

 glutTimerFunc(GAP, timer, minutes);
 glutPostRedisplay();

}


void displayClock()
{


   stringstream okk;
   okk<<hours;
   string str,str1;
   okk>>str;
   str = str + " seconds left";
   //set_mat_prop(1.0,0.0,0.0);
   drawString(boat_x+3,boat_y+10,-8.5,str);



}




void curveboat(void )
{
    //nt = 6
    //boat
    glPushMatrix();
    //glColor3f(1.0,1.0,1.0);
    set_mat_prop(0.0,0.0,0.0);
    glTranslatef(-0.5,-0.9,1.55);
    glRotatef(90,-1,0,0);
    glScalef(0.63,0.53,0.57);
    //glScalef(0.27,0.055,0.8);
    drawpyramid();
    glPopMatrix();


//right
    glPushMatrix();
    glTranslatef(0.8,-0.5,2);
    set_mat_prop(0.647, 0.165, 0.165);

    glRotatef( 90+90, 0.0, 0.0, 1.0);
    //glScalef(1.3,1.5,1.5);
    glScalef(0.7,0.5,1.1);
    //set_mat_prop(0.184, 0.310, 0.310);
    //set_mat_prop(0.0, 0.0, 0.0);
    bottleBezier(6);
    glPopMatrix();

   //left
    glPushMatrix();
    glTranslatef(-0.7,-0.5,2);

    //set_mat_prop(0.184, 0.310, 0.310);
    set_mat_prop(0.647, 0.165, 0.165);
    glScalef(0.7,0.5,1.1);
    bottleBezier(6);
    glPopMatrix();

    glPushMatrix();
    //set_mat_prop(0.961, 0.961, 0.961);
    set_mat_prop(0.0,0.0,0.0);
    //glTranslatef(0.1,0.1,0.9);
    glTranslatef(0.1,-0.5,0.5);
    glRotatef(-90,0,1,0);
    //glScalef(0.8,0.18,0.5);
    glScalef(0.8,0.3,0.6);
    //glutSolidCone(1,0.7,50,10);
    bottleBezier(10);
    glPopMatrix();

    glPushMatrix();
    set_mat_prop(0.827, 0.827, 0.827);
    glTranslatef(0,-0.5,3.1);
    glScalef(0.78,0.3,0.3);
    //set_mat_prop(0.184, 0.310, 0.310);
    glutSolidCube(1);
    glPopMatrix();

    glPushMatrix();
    //glColor3f(1.0,0.0,0.0);
    set_mat_prop(0.627, 0.322, 0.176);
    glTranslatef(0.04,-0.7,1.9);
    glScalef(0.27,0.055,0.8);
    glutSolidCube(3);
    glPopMatrix();



}
void boat(void)
{
    //boat

    set_mat_prop(0.0,0.0,0.0);
    glPushMatrix();

   glTranslatef(0,-2,0);
   //glRotatef(boat_rotate, 0,1,0);
    glScalef(0.51,0.1,1.5);
    cube(0.0,0.0,0.0);
    //glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    //pyramid drawing
    glPushMatrix();
    glTranslatef(0,-2,0);
    glRotatef(-50, 1,0,0);
    glScalef(0.26,0.30,0.1);

    drawpyramid();
    glPopMatrix();

    glPushMatrix();
   glTranslatef(0-0.01,-2,0+1.3);
   glRotatef(50, 1,0,0);
    glScalef(0.26,0.30,0.1);
    set_mat_prop(0.0,0.0,0.0);
    drawpyramid();
    glPopMatrix();
}

void flower1(void)
    {glPushMatrix();
    glTranslatef(-2,-0.28,14);
    glRotatef(60,0,1,0);
    glScalef(0,1,1);
    set_mat_prop(1.0,0.0,0.0);
    glutSolidSphere(0.5,5,2);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-2,-0.3,14);
    glRotatef(90,1,0,0);
    draw_cylinder(0.05, 0.5, 0.333, 0.420, 0.184);
    glPopMatrix();
    }

    void bird(void)
    {
        glPushMatrix();
    glTranslatef(3.0,1.5,11);
    glScalef(0.9,0.5,1.0);
    //cube(0.627, 0.322, 0.176);
    set_mat_prop(1.0, 1.0, 1.000);
    glutSolidIcosahedron();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.9,1.5,11);
    glScalef(0.9,0.5,1.0);
    //cube(0.627, 0.322, 0.176);
    set_mat_prop(1.0, 1.0, 1.000);
    glutSolidIcosahedron();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.5,1.5,11);
    glScalef(0.3,1.1,1.0);
    //cube(0.627, 0.322, 0.176);
    set_mat_prop(1.0, 1.0, 1.000);
    glutSolidIcosahedron();
    glPopMatrix();
    }

//cat drawing

void
makeLegAndPaw ()
{

	// Leg

	glPushMatrix ();
	glTranslatef (0, -0.7, 0);
	glScalef (0.1, 0.7, 0.1);
	glutSolidSphere (1, 50, 50);
	glPopMatrix ();

	// Paw

	glPushMatrix ();
	glTranslatef (0.05, -1.3, 0);
	glScalef (0.12, 0.1, 0.09);
	glutSolidSphere (1, 50, 50);
	glPopMatrix ();
}
void
clipTailBody ()
{
	double clip_plane0[] = { 0, -1, 0, 0 };
	double clip_plane1[] = { -1, 0, 0, 0 };
	glClipPlane (GL_CLIP_PLANE0, clip_plane0);
	glClipPlane (GL_CLIP_PLANE1, clip_plane1);

	glEnable (GL_CLIP_PLANE0);
	glEnable (GL_CLIP_PLANE1);
	glutSolidTorus (0.08, 0.4, 50, 50);
	glDisable (GL_CLIP_PLANE0);
	glDisable (GL_CLIP_PLANE1);
}

void
clipTailBall ()
{
	double clip_plane0[] = { -1, 0, 0, 0 };
	glClipPlane (GL_CLIP_PLANE0, clip_plane0);

	glEnable (GL_CLIP_PLANE0);

	glutSolidSphere (0.08, 50, 50);
	glDisable (GL_CLIP_PLANE0);
}

void
clipSmile ()
{

	double clip_plane0[] = { 0, -1, 0, 0 };
	glClipPlane (GL_CLIP_PLANE0, clip_plane0);

	glEnable (GL_CLIP_PLANE0);

	glutSolidTorus (0.01, 0.05, 50, 50);
	glDisable (GL_CLIP_PLANE0);
}

double front_legs_angle_param = -15;
double back_legs_angle_param	= 0;



void
generateCatto ()
{

	// Body

	glPushMatrix ();
	glTranslatef (0, 1.6, 0);
	glScalef (0.8, 0.6, 0.4);
	glutSolidSphere (1, 50, 50);
	glPopMatrix ();

	// Back legs

	glPushMatrix ();
	glTranslatef (-0.4, 1.7, 0.3);
	glRotatef (back_legs_angle_param, 0, 0, 1);
	makeLegAndPaw ();
	glPopMatrix ();

	glPushMatrix ();
	glTranslatef (-0.4, 1.7, -0.3);
	glRotatef (back_legs_angle_param, 0, 0, 1);
	makeLegAndPaw ();
	glPopMatrix ();

	// Front legs

	glPushMatrix ();
	glTranslatef (0.4, 1.7, 0.3);
	glRotatef (-front_legs_angle_param, 0, 0, 1);
	makeLegAndPaw ();
	glPopMatrix ();

	glPushMatrix ();
	glTranslatef (0.4, 1.7, -0.3);
	glRotatef (-front_legs_angle_param, 0, 0, 1);
	makeLegAndPaw ();
	glPopMatrix ();

	// Head
	glPushMatrix ();
	glTranslatef (0.8, 2, 0);
	glScalef (0.25, 0.4, 0.3);
	glutSolidSphere (1, 50, 50);
	glPopMatrix ();

	// Ears
	glPushMatrix ();
	glTranslatef (0.8, 2.3, -0.15);
	glRotatef (-30, 1, 0, 0);
	glScalef (0.3, 1.3, 1);
	glRotatef (90, 0, 0, 1);
	glScalef (0.2, 0.7, 0.2);
	glutSolidTetrahedron ();
	glPopMatrix ();

	glPushMatrix ();
	glTranslatef (0.8, 2.3, 0.15);
	glRotatef (30, 1, 0, 0);
	glScalef (0.3, 1.3, 1);
	glRotatef (90, 0, 0, 1);
	glScalef (0.2, 0.7, 0.2);
	glutSolidTetrahedron ();
	glPopMatrix ();

	// Tail
	glPushMatrix ();
	glTranslatef (-0.7, 2.2, 0);
	glRotatef (0, 0, 0, 1);
	clipTailBody ();
	glPopMatrix ();

	glPushMatrix ();
	glTranslatef (-1.5, 2.2, 0);
	glRotatef (180, 0, 0, 1);
	clipTailBody ();
	glPopMatrix ();

	glPushMatrix ();
	glTranslatef (-1.5, 2.6, 0);
	clipTailBall ();
	glPopMatrix ();


	// Smile
	set_mat_prop (0, 0, 0);

	glPushMatrix ();
	glTranslatef (1.03, 1.9, 0.05);
	glRotatef (70, 0, 1, 0);
	clipSmile ();
	glPopMatrix ();

	glPushMatrix ();
	glTranslatef (1.03, 1.9, -0.05);
	glRotatef (110, 0, 1, 0);
	clipSmile ();
	glPopMatrix ();

	// Eyes

	glPushMatrix ();
	glTranslatef (0.99, 2.2, -0.08);
	glRotatef (30, 0, 0, 1);
	glRotatef (100, 0, 1, 0);
	glScalef (1.3, 1, 1);
	clipSmile ();
	glPopMatrix ();

	glPushMatrix ();
	glTranslatef (0.99, 2.2, 0.08);
	glRotatef (30, 0, 0, 1);
	glRotatef (80, 0, 1, 0);
	glScalef (1.3, 1, 1);
	clipSmile ();
	glPopMatrix ();

	// Eyeballs

	glPushMatrix ();
	glTranslatef (1.025, 2.115, 0.095);
	glScalef (0.01, 0.035, 0.01);
	glutSolidSphere (1, 50, 50);
	glPopMatrix ();

	glPushMatrix ();
	glTranslatef (1.025, 2.115, -0.095);
	glScalef (0.01, 0.035, 0.01);
	glutSolidSphere (1, 50, 50);
	glPopMatrix ();
}

void tree1()
{
    //tree

    glPushMatrix();
    glTranslatef(3.5,1,11);
    glRotatef(90,1,0,0);
    //glScalef(1,2,1);
    draw_cylinder(0.16, 5.5, 0.627, 0.322, 0.176);
    glPopMatrix();

    //tree leaf
    glPushMatrix();
    glTranslatef(3.0,1.5,11);
    glScalef(0.9,0.5,1.0);
    //cube(0.627, 0.322, 0.176);
    set_mat_prop(0.0, 1.0, 0.000);
    glutSolidIcosahedron();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.9,1.5,11);
    glScalef(0.9,0.5,1.0);
    //cube(0.627, 0.322, 0.176);
    set_mat_prop(0.0, 1.0, 0.000);
    glutSolidIcosahedron();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.5,2,11);
    glScalef(0.4,0.8,1.0);
    //cube(0.627, 0.322, 0.176);
    set_mat_prop(0.0, 1.0, 0.000);
    glutSolidDodecahedron();
    glPopMatrix();



}

void scene1( void )
{
    GLfloat eyeX = 2.5;
    GLfloat eyeY = 2.5;
    GLfloat eyeZ = 20;

    GLfloat lookX = 1;
    GLfloat lookY = -0.5;
    GLfloat lookZ = -11;
//displayClock();
     //timer(0);
    glPushMatrix();
    glTranslatef(-1,4,0);
    bird();
    glPopMatrix();


    glPushMatrix();
    glTranslatef(-10,0,0);
    tree1();
    glPopMatrix();
    //glFlush();
    string str1,str2,str3,str4,str5,str6,str7,str8,str9,str10;
    glPushMatrix();
    str1= "Adventure in the River of Sundarbans";
    glScalef(5.0,5.0,3.0);
    set_mat_prop(0.0,0.0,0.0);
    drawString(-1.5,1,0,str1);
    glPopMatrix();
    glPushMatrix();
    set_mat_prop(0.500, 0.0, 0.000);
    glScalef(2.0,2.0,2.0);
    str3= "Key Instructions:";
    drawString(3,0,0,str3);
    set_mat_prop(0.294, 0.000, 0.510);
    str2= "New Game  v" ;
    drawString(3.5,-2,0,str2);

    str4= "Forward      Up-Arrow " ;
    drawString(3.5,-3,0,str4);
    str5= "Back      Down-Arrow" ;
    drawString(3.5,-4,0,str5);
    str6= "Left      Left-Arrow" ;
    drawString(3.5,-5,0,str6);
    str7= "Right      Right-Arrow" ;
    drawString(3.5,-6,0,str7);
    str8= "Rotate-Right      4" ;
    drawString(3.5,-7,0,str8);
    str9= "Rotate-Left      5" ;
    drawString(3.5,-8,0,str9);
    str10= "Exit      0" ;
    drawString(3.5,-9,0,str10);
    glPopMatrix();


}


//tree

#define M_PI		3.14159265358979323846
void lineTo(float x, float y, float len, float angle)
{
    glBegin(GL_LINES);
    glVertex2f(x, y);
    glVertex2f(x+len*cos(angle), y+len*sin(angle));
    glEnd();
}

void tree(float x, float y, float len, float angle, float len_div, float angle_dif, int depth)
{
    lineTo(x, y, len, angle);
    if(depth == 0) return;

    // left
    tree(x+len*cos(angle), y+len*sin(angle), len*len_div, angle-angle_dif, len_div, angle_dif, depth-1);

    // right
    tree(x+len*cos(angle), y+len*sin(angle), len*len_div, angle+angle_dif, len_div, angle_dif, depth-1);
}
//tree4
int flag1 = 1;
GLfloat angle = 0.0f;
GLfloat angle2 = 0.0f;
int moving, startx, starty;
float r = 0.0f;//0.47f;
float g = 1.0f;//0.0f;
float b = 0.0f;//0.74f;
void makecylinder(float height,float Base)
{
	GLUquadricObj *qobj;
	qobj = gluNewQuadric();
	glColor3f(r, g, b);
	glPushMatrix();
	glRotatef(-90, 1.0f, 0.0f, 0.0f);
	gluCylinder(qobj, Base, Base - (0.2 * Base), height, 20, 20);
	glPopMatrix();
}
void maketree(float height,float Base)
{

	glPushMatrix();
	float angle;
	makecylinder(height, Base);
	glTranslatef(0.0f, height,0.0f);
	height -= height * 0.2f;
	Base -= Base * 0.3f;

	if (height > 1)
	{
		angle = 22.5f;
		glPushMatrix();
		glRotatef(angle, -1.0f, 0.0f, 0.0f);
		maketree(height, Base);
		glPopMatrix();
		glPushMatrix();
		glRotatef(angle, 0.5f, 0.0f, 0.866f);
		maketree(height, Base);
		glPopMatrix();
		glPushMatrix();
		glRotatef(angle, 0.5f, 0.0f, -0.866f);
		maketree(height, Base);
		glPopMatrix();
	}
	glPopMatrix();
}

void newflower()
{
    //insect
    glPushMatrix();
    glTranslatef(4,0.3,16.8);
    //glRotatef(90,0,0,1);
    glScalef(0.07,0.07,0.07);
    set_mat_prop(0.863, 0.078, 0.235);
    insectBezier(1);
    glPopMatrix();
}


void tree5(void)
{
  glPushMatrix();
    //glViewport(0,0,700,700);
    glTranslatef(7,-2,8);
    glScalef(0.7,0.7,0.7);
    set_mat_prop(0.502, 0.000, 0.000);
    maketree(2,0.2);
    glPopMatrix();
    //leaf1
    glPushMatrix();
    glTranslatef(6,2,8);
    glScalef(1,0.7,0.7);
    set_mat_prop(0.0, 0.539, 0.000);
    glutSolidSphere(0.5,5,5);
    glPopMatrix();
    //leaf2
    glPushMatrix();
    glTranslatef(5.7,1.7,8);
    glScalef(1,1,0.7);
    set_mat_prop(0.0, 0.539, 0.000);
    glutSolidSphere(0.5,10,8);
    glPopMatrix();
    //leaf3
    glPushMatrix();
    glTranslatef(5.5,1,8.5);
    glRotatef(-90,0,0,1);
    glScalef(0.2,0.3,0.55);
    set_mat_prop(0.0, 0.539, 0.000);
    glutSolidDodecahedron();
    glPopMatrix();
    //leaf4
     glPushMatrix();
    glTranslatef(6.7,1.8,8);
    glRotatef(-40,0,0,1);
    glScalef(0.2,0.3,0.5);
    set_mat_prop(0.0, 0.539, 0.000);
    glutSolidDodecahedron();
    glPopMatrix();
    //leaf5
    glPushMatrix();
    glTranslatef(7.3,1.8,9);
    glRotatef(-50,0,0,1);
    glScalef(0.2,0.3,0.5);
    set_mat_prop(0.0, 0.539, 0.000);
    glutSolidDodecahedron();
    glPopMatrix();
}
void collidingRocks(void)
{

    //rock7
    glPushMatrix();
    glTranslatef(23,-3,5);
    set_mat_prop(0.412, 0.412, 0.412);
    glutSolidDodecahedron();
    glPopMatrix();

     //rock6
    glPushMatrix();
    glTranslatef(20,-2,-5);
    set_mat_prop(0.412, 0.412, 0.412);
    glutSolidDodecahedron();
    glPopMatrix();
    //rock5
    glPushMatrix();
    glTranslatef(12,-2,3);
    glScalef(0.7,0.7,0.7);
    set_mat_prop(0.412, 0.412, 0.412);
    glutSolidDodecahedron();
    glPopMatrix();
    //rock4
    glPushMatrix();
    glTranslatef(3,-2,-3);
    set_mat_prop(0.412, 0.412, 0.412);
    glutSolidDodecahedron();
    glPopMatrix();

    //rock3
    glPushMatrix();
    glTranslatef(4,-2,3);
    set_mat_prop(0.412, 0.412, 0.412);
    glutSolidDodecahedron();
    glPopMatrix();
    //rock2
    glPushMatrix();
    glTranslatef(-2.5,-2,6);
    set_mat_prop(0.412, 0.412, 0.412);
    glutSolidDodecahedron();
    glPopMatrix();
    //rock1
    glPushMatrix();
    glTranslatef(-2.8,-2,16);
    set_mat_prop(0.412, 0.412, 0.412);
    glutSolidDodecahedron();
    glPopMatrix();
    //boat
    glPushMatrix();
    glTranslatef(10,-1,-10);
    glRotatef(-80,0,1,0);
    glScalef(1.2,1.2,1.2);
    //set_mat_prop(1.0, 1.0, 0.0);
    //glutSolidDodecahedron();
    glColor3f(1.0,0.0,0.0);
    boat();
    glPopMatrix();
}
int wincheck=0;
int colicheck=0;
int collision(void)
{
    //left sky+grass
    if((boat_x <=-5) && (boat_z>= -11 && boat_z<= 20))
        colicheck =1;
    //middle sky
    if(boat_z<=-10.5)
        colicheck =1;
    //right grass
    if((boat_x >=2.8 && boat_x<=10) && (boat_z>= 6 && boat_z<= 18))
        colicheck =1;
    //right sky
    if((boat_x >22.5) && (boat_z>= -11 && boat_z<= 20))
        colicheck =1;

    //boat
    if((boat_x >= 8 && boat_x<= 15) && (boat_z>= -10.5 && boat_z<= -8))
        colicheck =1;
    //1st rock
     if((boat_x >= -3.5 && boat_x<= 0) && (boat_z>=14.5 && boat_z<=17.5))
        colicheck =1;
    //2nd rock
    if((boat_x >= -3.5 && boat_x<=0.3) && (boat_z>=4.5 && boat_z<=7))
        colicheck =1;
    //3rd rock
    if((boat_x >=2 && boat_x<=4.5) && (boat_z>=2 && boat_z<=4))
        colicheck =1;
    //4th rock
    if((boat_x >=1.5 && boat_x<=3.5) && (boat_z>=-4 && boat_z<=-1))
        colicheck =1;
    //5th rock
    if((boat_x >=10.5 && boat_x<=13.5) && (boat_z>=1.5 && boat_z<=4.5))
        colicheck =1;
    //6th rock
    if((boat_x >=17 && boat_x<=21.5) && (boat_z>=-8 && boat_z<=-4))
        colicheck =1;
    //7th rock
    if((boat_x >=20 && boat_x<=23.5) && (boat_z>=3 && boat_z<=6))
        colicheck =1;
    //Finishing line
    if((boat_x >=12 && boat_x<=23) && (boat_z>=12.5 && boat_z<=14.5))
        { wincheck=1;
        colicheck =1;}


    if(colicheck==1)
        {
            return 1;
            sceneCall=2;
         }
    else
        return 0;

}

void scene3( void )
{
        eyeX = 2.5;
        eyeY = 2.5;
        eyeZ = 20;

        lookX = 1;
        lookY = -0.5;
        lookZ = -11;

        boat_x=2;
        boat_y= -2.5;
        boat_z=14.5;

    string str1,str2,str3,str4;

    if(wincheck==0)
    {
    glPushMatrix();
    str1= "Game Over";
    glScalef(1.0,1.0,1.0);
    set_mat_prop(1.0,0.0,0.0);
    drawString(-3,0,0,str1);
    glPopMatrix();

    glPushMatrix();
    str2="New Game v";
    glScalef(1.0,1.0,1.0);
    set_mat_prop(0.502, 0.000, 0.502);
    drawString(-3,-3,0,str2);
    glPopMatrix();

    glPushMatrix();
    str3="Home Page  n";
    glScalef(1.0,1.0,1.0);
    set_mat_prop(0.502, 0.000, 0.502);
    drawString(-3,-4.5,0,str3);
    glPopMatrix();

    glPushMatrix();
    str4="Exit  0";
    glScalef(1.0,1.0,1.0);
    set_mat_prop(0.502, 0.000, 0.502);
    drawString(-3,-6,0,str4);
    glPopMatrix();


    }
    else if(wincheck==1)
    {
    glPushMatrix();
    str1= "You Won!";
    glScalef(1.0,1.0,1.0);
    set_mat_prop(1.0,0.0,0.0);
    drawString(-3,0,0,str1);
    glPopMatrix();

    glPushMatrix();
    str2="New Game v";
    glScalef(1.0,1.0,1.0);
    set_mat_prop(0.502, 0.000, 0.502);
    drawString(-3,-3,0,str2);
    glPopMatrix();

    glPushMatrix();
    str3="Home Page  n";
    glScalef(1.0,1.0,1.0);
    set_mat_prop(0.502, 0.000, 0.502);
    drawString(-3,-4.5,0,str3);
    glPopMatrix();

    glPushMatrix();
    str4="Exit  0";
    glScalef(1.0,1.0,1.0);
    set_mat_prop(0.502, 0.000, 0.502);
    drawString(-3,-6,0,str4);
    glPopMatrix();

    }


}

void scene2( void )
{
     //collision();
    glPushMatrix();
    collidingRocks();
    glPopMatrix();
    //left sky
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, ID[2]);
    glTranslatef(-15,-5,-11);
    //glRotatef(90,0,1,0);
    glScalef(1,40,28);
    cube();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();


//right sky

glPushMatrix();

   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, ID[2]);
    glTranslatef(24.2,-5,-11);
    glScalef(1,40,28);
    cube();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();


 // Middle sky
   glPushMatrix();

   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, ID[0]);
    glTranslatef(-17,-5,-11);
    glScalef(48,40,0.05);
    cube();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

//left grass
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, ID[2]);
   glTranslatef(-14,-4,-10);
    glScalef(7,0.5,27);
    cube();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

//right water
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, ID[1]);
   glTranslatef(-8,-4,-10.9);
    glScalef(38,0.1,27.9);
    cube();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
//boat
    glPushMatrix();
    glTranslatef(boat_x,boat_y,boat_z);

    //glRotatef(180,0,1,0);
    //glRotatef(-72, 0,0,1);
    glRotatef(boat_rotate, 0,1,0);
    glScalef(0.85,0.85,0.85);
    curveboat();

    glPopMatrix();



    //timer
    glPushMatrix();
    //set_mat_prop(0.0,0.0,1.0);
    timer(0);
    glFlush();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(4,-4,10);
    glScalef(5,2,6);
    set_mat_prop(0.000, 0.392, 0.000);
    bottleBezier(10);
    glPopMatrix();

/*
    glPushMatrix();
    flower1();
    glPopMatrix();
    glPushMatrix();
    glTranslatef(-0.2,0,0.4);
    flower1();
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0.4,0,0.8);
    flower1();
    glPopMatrix();
*/

    //bird1
    glPushMatrix();
    glScalef(0.5,0.5,0.5);
    glTranslatef(-2,18,-7);
    bird();
    //bird2
    glPopMatrix();
    glPushMatrix();
    glScalef(0.5,0.5,0.5);
    glTranslatef(-2.5,20,-7);
    bird();
    glPopMatrix();

   /*
    //brezier flower
    glPushMatrix();
    //glRotatef(90,0,0,1);
    newflower();
    glPopMatrix();

    //flower2
    glPushMatrix();
    glTranslatef(0.72,0,-0.3);
    newflower();
    glPopMatrix();
    //flower3
    glPushMatrix();
    glTranslatef(0.5,0,-0.6);
    newflower();
    glPopMatrix();
   */

    //tree with decohydron and cylinder
    glPushMatrix();
    glTranslatef(20,0,-10);
    tree1();
    glPopMatrix();
    //2d fractal trees
    glPushMatrix();
    glTranslatef(10,0,13);
    glRotatef(180,1,0,0);
    glScalef(3,3,3);
    set_mat_prop(0.333, 0.420, 0.184);
    tree(3, 2, 1, -M_PI/2.0, 0.8, M_PI/12.0, 15);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(7,0,15.5);
    glRotatef(180,1,0,0);
    glScalef(3,3,3);
    set_mat_prop(0.333, 0.420, 0.184);
    tree(3, 2, 1, -M_PI/2.0, 0.8, M_PI/12.0, 15);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,0,16);
    glRotatef(180,1,0,0);
    glScalef(3,3,3);
    set_mat_prop(0.333, 0.420, 0.184);
    tree(3, 2, 1, -M_PI/2.0, 0.8, M_PI/12.0, 15);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,0,14);
    glRotatef(180,1,0,0);
    glScalef(3,3,3);
    set_mat_prop(0.333, 0.420, 0.184);
    tree(3, 2, 1, -M_PI/2.0, 0.8, M_PI/12.0, 15);
    glPopMatrix();


    //Finishing Line
    glPushMatrix();
    glTranslatef(12,0,12);
    glScalef(20,0.7,0.3);
    cube(0.000, 0.392, 0.000);
    glPopMatrix();

    //3d fractal trees
    glPushMatrix();
    tree5();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.5,0,3.5);
    glRotatef(30,0,-1,0);
    tree5();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,-3.8,11);
    glScalef(0.7,0.5,0.5);
    set_mat_prop(0.000, 0.392, 0.000);
    maketree(2,0.17);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.5,-3.8,13);
    glScalef(0.8,0.5,0.5);
    set_mat_prop(0.000, 0.392, 0.000);
    maketree(2,0.17);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.3,-3.8,14.5);
    glScalef(0.8,0.5,0.5);
    set_mat_prop(0.000, 0.392, 0.000);
    maketree(2,0.17);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-8,-3,11.5);
    glScalef(1,1,1);
    set_mat_prop(0.000, 0.392, 0.000);
    maketree(4,0.17);
    glPopMatrix();
    //2d curve grass flowers
    glPushMatrix();
    glTranslatef(7,-2,11);
    glRotatef(180,1,0,0);
    glScalef(0.7,0.7,0.7);
    set_mat_prop(0.502, 0.502, 0.000);
    tree(3, 2, 1, -M_PI/2.0, 0.8, M_PI/12.0, 15);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,-2,9);
    glRotatef(180,1,0,0);
    glScalef(0.7,0.7,0.7);
    set_mat_prop(0.502, 0.000, 0.502);
    tree(3, 2, 1, -M_PI/2.0, 0.8, M_PI/12.0, 15);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(4.5,-2,10);
    glRotatef(180,1,0,0);
    glScalef(0.7,0.7,0.7);
    set_mat_prop(0.604, 0.804, 0.196);
    tree(3, 2, 1, -M_PI/2.0, 0.8, M_PI/12.0, 15);
    glPopMatrix();



    glPushMatrix();
    glTranslatef(-8,-3,5);
    glScalef(0.8,0.8,0.8);
    set_mat_prop(0.741, 0.718, 0.420);
    glRotatef(rotation_test,0,0.3,0);
    generateCatto();
    glPopMatrix();

}

static void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glClearColor(1, 1, 1,1);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-4, 4, -3, 3, 3.0, 200.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eyeX,eyeY,eyeZ, lookX,lookY,lookZ, 0,1,0);

    glRotatef(rot, 0,1,0);


int flagcoli=0;


if ( sceneCall == 0)
        scene1();
else if ( sceneCall == 1 )
{
 scene2();
 flagcoli = collision();
 if(flagcoli==0)
    sceneCall=1;
else
    sceneCall=2;
}
else if ( sceneCall == 2 )
        {scene3();
        }

//    // Axes
    //axes();

    glFlush();
    glutSwapBuffers();
}
void boat_init(void)
{
        eyeX = 2.5;
        eyeY = 2.5;
        eyeZ = 20;

        lookX = 1;
        lookY = -0.5;
        lookZ = -11;

        boat_x=2;
        boat_y= -2.5;
        boat_z=14.5;
        boat_rotate=0;
        wincheck=0;
        colicheck=0;
        hours=120;
        j=1;
        sceneCall =1;
}

static void key(unsigned char key, int x, int y)
{
    float x_, z_, r, theta, dx, dz, dx_norm, dz_norm, r_=1, turn_angle_step=5, height_diff_one_less, height_diff_thresh_dist;

    x_=lookX-eyeX;
    z_=lookZ-eyeZ;
    r=sqrt(x_*x_+z_*z_);

    if(x_==0)
        theta = 90;
    else
        theta=atan(z_/x_) * 180 / PI;

    if((z_>0 && theta<0) || (z_<0 && theta>0))
        theta += 180;

    dx = r_*cos(theta * PI / 180);
    dz = r_*sin(theta * PI / 180);

    dx_norm = r_*cos((theta-90) * PI / 180);
    dz_norm = r_*sin((theta-90) * PI / 180);

    switch (key)
    {

    case '0':
        exit(0);
        break;

    // moving the look at point at a circular path or up-down
    case 'a':
        theta-=turn_angle_step;
        theta = theta * PI / 180;

        lookX=r*cos(theta)+eyeX;
        lookZ=r*sin(theta)+eyeZ;
        break;
    case 'w':
        lookY++;
        break;
    case 's':
        lookY--;
        break;
    case 'd':
        theta+=turn_angle_step;
        theta = theta * PI / 180;

        lookX=r*cos(theta)+eyeX;
        lookZ=r*sin(theta)+eyeZ;
        break;

    // Moving the camera front-back-left-right
    case 'j':
            eyeX += dx_norm;
            eyeZ += dz_norm;

            lookX += dx_norm;
            lookZ += dz_norm;
        break;
    case 'i':
            eyeX += dx;
            eyeZ += dz;

            lookX += dx;
            lookZ += dz;


        break;
    case 'k':
            eyeX -= dx;
            eyeZ -= dz;

            lookX -= dx;
            lookZ -= dz;

        break;
    case 'l':
            eyeX -= dx_norm;
            eyeZ -= dz_norm;

            lookX -= dx_norm;
            lookZ -= dz_norm;
        break;

    // Moving the camera up-down
    case '+':
        eyeY++;
        lookY++;
        break;
    case '-':
        eyeY--;
        lookY--;
        break;


    // rotating the whole scene
    case ',':
        rot--;
        break;
    case '.':
        rot++;
        break;

    // rotating the camera around the look at point
    case 'g':
        theta += 180;
        theta += turn_angle_step;
        theta = theta * PI / 180;

        eyeX=r*cos(theta)+lookX;
        eyeZ=r*sin(theta)+lookZ;
        break;
    case 'h':
        theta += 180;
        theta -= turn_angle_step;
        theta = theta * PI / 180;

        eyeX=r*cos(theta)+lookX;
        eyeZ=r*sin(theta)+lookZ;
        break;

    // look far or near
    case 'r':
        lookX += dx;
        lookZ += dz;
        break;
    case 'f':
        lookX -= dx;
        lookZ -= dz;
        break;

    // on-off switches to whichever it is attached
    case '1':
    {
        right_light = !right_light;
        light_right(20,5,-5);

        break;

    }
    case '2':
    {
        left_light = !left_light;
        light_left(-20,10,0);

        break;

    }
    case '3':
    {
        spot_light = !spot_light;
        light_spot(0,10,10);

        break;

    }
    case 'z':

        activate_amb = !activate_amb;
        light_right(20,5,-5);
        light_left(20,10,0);
        light_spot(0,10,10);;
        break;


    case 'x':

        activate_spec = !activate_spec;
        light_right(20,5,-5);
        light_left(20,10,0);
        light_spot(0,10,10);;
        break;

    case 'c':

        activate_diff = !activate_diff;
        light_right(20,5,-5);
        light_left(20,10,0);
        light_spot(0,10,10);;
        break;



      case '9':
        cat_rotate=true;
        //rotation_test=0;
        break;


    case '4':
        boat_rotate -= 1;
        break;
    case '5':
        boat_rotate += 1;
        break;


    case 'v':
        //flagcoli=0;
        boat_init();

        if(j==0)
            j=1;
        else
            j=0;
        break;
    case 'n':
        sceneCall = 0;
        break;



    }

    glutPostRedisplay();
}

void animate()
{
    if (cat_rotate == true)
    {
        rotation_test= rotation_test-1;
        if(rotation_test<= -90)
            rotation_test=0;

    }
    glutPostRedisplay();
}


void specialkey (int key, int x, int y)
{
    switch(key)
    {
    case GLUT_KEY_UP:
        boat_z -= 0.2;
        lookZ -= 0.1;
        eyeZ -= 0.1;
        //collision();
        glutPostRedisplay();
        break;

    case GLUT_KEY_DOWN:
        boat_z += 0.2;
        lookZ += 0.02;
        eyeZ += 0.02;
        //collision();
        glutPostRedisplay();
        break;
    case GLUT_KEY_LEFT:
        boat_x -= 0.1;
        //collision();
        glutPostRedisplay();
        break;
    case GLUT_KEY_RIGHT:
        boat_x += 0.1;
        //collision();
        glutPostRedisplay();
        break;


    }
}


int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowSize(width,height);
    glutInitWindowPosition(30,30);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("1607009 project");

    glutDisplayFunc(display);
    glutKeyboardFunc(key);

    glutReshapeFunc(resize);
    glutSpecialFunc(specialkey);

    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING);


    LoadTexture("I:\\Regular Class\\4-2\\Graphics Lab\\Texture\\jungle2.bmp",0);
    LoadTexture("I:\\Regular Class\\4-2\\Graphics Lab\\Texture\\water2.bmp",1);
    LoadTexture("I:\\Regular Class\\4-2\\Graphics Lab\\Texture\\jungle4.bmp",2);
    LoadTexture("I:\\Regular Class\\4-2\\Graphics Lab\\Texture\\grass.bmp",5);



    printf("Use 'w' to look up,\n 's' to look down,\n 'd' to look right,\n and 'a' to look left.\n\n");
    printf("Use 'i' to move camera front (i.e. zoom in),\n 'k' to move camera back (i.e. zoom out),\n 'l' to move camera right,\n and 'j' to move camera left.\n\n");
    printf("Use '+' to move camera up\n and '-' to move camera down.\n\n");
    printf("Use 'g' to rotate left,\n and 'h' to rotate right taking the look at point as center.\n\n");
    printf("Use 'r' to look far,\n and 'f' to look near.\n\n");
    printf("Use '1' and '2' to toggle tube lights switches and 3 for spotlight\n\n");
    printf("Use 'z' to toggle ambient,\n 'c' to toggle diffusion,\n and 'x' to toggle specular light property for all lights.\n\n\n\n");
    glutIdleFunc(animate);
    glClearColor(0.690, 0.769, 0.871,1);
    //collision();
    glutMouseFunc(processMouse);
    glutMainLoop();


    return EXIT_SUCCESS;
}
