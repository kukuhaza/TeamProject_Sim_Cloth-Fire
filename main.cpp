
#define STB_IMAGE_IMPLEMENTATION

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>

#include "stb_image.h";
#include <glut.h>
#include <GL\GLU.h>
#include <glm\glm.hpp>


using namespace glm;

#define NUMBER_OF_FIRE 1000


float timeStep;

//key control variable
int rotate[3] = { 0 };
bool m_start = true;
bool isrender = false;
bool master = true;
bool strip = false;

const float LENGTH = 60.0f;
vec3 gravity_point = vec3(0.0f, -9.8f, 0.0f);
vec3 ground = vec3(0.0f, -10 - (LENGTH / 2), 0.0f);


struct particle {
	vec3 Position, Velocity;
	vec3 Color;

	bool active;
	float Life;
	float mass;
	float fade;
	//half length of side
	float hl;

	void init() {
	
		mass = 0.01f;
		active = true;
		
		//randomly create particle position in top of cube
		//Position = vec3(rand() % 21 - 10, rand()%21-10-(LENGTH/2) , rand() % 20 - 10);
		Position = vec3(rand() % 21 - 10, 10 - (LENGTH / 2), rand() % 20 - 10);
		Life = 1.0f;
		//randomly set fade value -> minimum=0.007f
		fade = ((rand() % 100) / 1000.0f) + 0.007f;
		hl = 4.0f;

		//shoot to upper side -> velocity. y > 0
		Velocity.x = rand() % 50 - 25;
		Velocity.y = rand() % 50;
		Velocity.z = rand() % 50 - 25;
	}

	void Movement(float t, vec3 force) {
		vec3 acc = force / mass;
		Velocity = Velocity + acc * t;
		Position = Position + Velocity * t;
	}
	
};


GLuint mTexture[2];
bool LoadMeshFromFile(const char* filename, int k)
{

	glGenTextures(1, &mTexture[k]);

	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		printf("ERROR: NO File. failed to bind texture\n");
		return false;
	}
	int width, height, channel;
	unsigned char *image = stbi_load_from_file(fp, &width, &height, &channel, 4);
	fclose(fp);

	glBindTexture(GL_TEXTURE_2D, mTexture[k]);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);

	return true;
}

class particle_system {
public:
	std::vector<particle> particles;

	void init(int n) {
		particles.clear();
		for (int i = 0; i < n; i++)
		{
			particle temp;    //create particle
			particles.push_back(temp);   //push in particles
		}
		for (int i = 0; i < particles.size(); i++) {
			particles[i].init();
		}
	}

	void Movement(float time) {
		for (int i = 0; i < particles.size(); i++) {
			//move particles to y axis faster
			vec3 force = vec3(0.0f, 0.8f, 0.0f);
			particles[i].Movement(time, force);
		}
	}


	void draw() {
		glEnable(GL_BLEND);
		glEnable(GL_TEXTURE_2D);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);

		glBindTexture(GL_TEXTURE_2D, mTexture[0]);

		for (int i = 0; i < particles.size(); i++) {
			if (particles[i].active) {

				vec3 pos = particles[i].Position;
				vec3 color = particles[i].Color;
				float hl = particles[i].hl;
				//use life for transparency
				glColor4f(color.x, color.y, color.z, particles[i].Life);

				//draw three Quads -> to see every angle
				glBegin(GL_QUADS);
				glTexCoord2f(0.0f, 0.0f);
				glVertex3f(pos.x - hl, pos.y + hl, pos.z);
				glTexCoord2f(0.0f, 1.0f);
				glVertex3f(pos.x - hl, pos.y - hl, pos.z);
				glTexCoord2f(1.0f, 1.0f);
				glVertex3f(pos.x + hl, pos.y - hl, pos.z);
				glTexCoord2f(1.0f, 0.0f);
				glVertex3f(pos.x + hl, pos.y + hl, pos.z);
				
				glTexCoord2f(0.0f, 0.0f);
				glVertex3f(pos.x , pos.y + hl, pos.z-hl);
				glTexCoord2f(0.0f, 1.0f);
				glVertex3f(pos.x, pos.y - hl, pos.z-hl);
				glTexCoord2f(1.0f, 1.0f);
				glVertex3f(pos.x,  pos.y - hl, pos.z+hl);
				glTexCoord2f(1.0f, 0.0f);
				glVertex3f(pos.x , pos.y + hl, pos.z+hl);

				glTexCoord2f(0.0f, 0.0f);
				glVertex3f(pos.x + hl, pos.y , pos.z - hl);
				glTexCoord2f(0.0f, 1.0f);
				glVertex3f(pos.x - hl, pos.y , pos.z - hl);
				glTexCoord2f(1.0f, 1.0f);
				glVertex3f(pos.x - hl, pos.y , pos.z + hl);
				glTexCoord2f(1.0f, 0.0f);
				glVertex3f(pos.x + hl, pos.y , pos.z + hl);
				glEnd();


			}
		}
		glDisable(GL_BLEND);
		glDisable(GL_TEXTURE_2D);
	}

	void Coloring() {
		for (int i = 0; i < particles.size(); i++) {
			//dead -> newly Init
			if (particles[i].Life <= 0.0f)
				particles[i].init();

			if (particles[i].Life < 0.6)
				particles[i].Color = vec3(1.0f, 0.0f, 0.0f); //red
			else if (particles[i].Life < 0.7)
				particles[i].Color = vec3(1.0f, 0.5f, 0.0f);//orange?
			else if (particles[i].Life < 0.8)
				particles[i].Color = vec3(1.0f, 1.0f, 0.0f);//yellow
			else if (particles[i].Life < 0.95)
				particles[i].Color = vec3(1.0f, 1.0f, 1.0f);//white

			//reduce Life
			particles[i].Life -= particles[i].fade;
			//gradually grow
			particles[i].hl += (particles[i].fade * 2);

		}
	}
};
particle_system ParticleSystem;


//class for cube
class rigid_cube {
public:
	vec3 Position;
	float hl;  //half length of side
	vec3 Color = vec3(1.0f, 1.0f, 1.0f);

	rigid_cube(vec3 pos, float num) {
		Position = pos;
		hl = num;
	}
	
	void draw() {

		glBegin(GL_QUADS);
		//bottom
		glColor3f(Color.x, Color.y, Color.z);
		glVertex3f(Position.x - hl, Position.y - hl, Position.z - hl);
		glVertex3f(Position.x - hl, Position.y - hl, Position.z + hl);
		glVertex3f(Position.x + hl, Position.y - hl, Position.z + hl);
		glVertex3f(Position.x + hl, Position.y - hl, Position.z - hl);

		//back
		
		glVertex3f(Position.x + hl, Position.y - hl, Position.z - hl);
		glVertex3f(Position.x - hl, Position.y - hl, Position.z - hl);
		glVertex3f(Position.x - hl, Position.y + hl, Position.z - hl);
		glVertex3f(Position.x + hl, Position.y + hl, Position.z - hl);


		//right
		glColor3f(0.0f, 1.0f, 1.0f);
		glVertex3f(Position.x + hl, Position.y + hl, Position.z - hl);
		glVertex3f(Position.x + hl, Position.y - hl, Position.z - hl);
		glVertex3f(Position.x + hl, Position.y - hl, Position.z + hl);
		glVertex3f(Position.x + hl, Position.y + hl, Position.z + hl);


		//left
		glVertex3f(Position.x - hl, Position.y + hl, Position.z - hl);
		glVertex3f(Position.x - hl, Position.y - hl, Position.z - hl);
		glVertex3f(Position.x - hl, Position.y - hl, Position.z + hl);
		glVertex3f(Position.x - hl, Position.y + hl, Position.z + hl);

		//front
		glColor3f(Color.x, Color.y, Color.z);
		glVertex3f(Position.x + hl, Position.y - hl, Position.z + hl);
		glVertex3f(Position.x - hl, Position.y - hl, Position.z + hl);
		glVertex3f(Position.x - hl, Position.y + hl, Position.z + hl);
		glVertex3f(Position.x + hl, Position.y + hl, Position.z + hl);
	
		glColor3f(1.0f, 0.0f, 0.0f);
		//top
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(Position.x - hl, Position.y + hl, Position.z - hl);
		glVertex3f(Position.x - hl, Position.y + hl, Position.z + hl);
		glVertex3f(Position.x + hl, Position.y + hl, Position.z + hl);
		glVertex3f(Position.x + hl, Position.y + hl, Position.z - hl);

		glEnd();

		//white line. activate with key 'L'
		if (strip) {
			glBegin(GL_LINE_STRIP);
			//bottom
			glColor3f(1.0f, 1.0f, 1.0f);
			glVertex3f(Position.x - hl, Position.y - hl, Position.z - hl);
			glVertex3f(Position.x - hl, Position.y - hl, Position.z + hl);
			glVertex3f(Position.x + hl, Position.y - hl, Position.z + hl);
			glVertex3f(Position.x + hl, Position.y - hl, Position.z - hl);

			//back
			glVertex3f(Position.x + hl, Position.y - hl, Position.z - hl);
			glVertex3f(Position.x - hl, Position.y - hl, Position.z - hl);
			glVertex3f(Position.x - hl, Position.y + hl, Position.z - hl);
			glVertex3f(Position.x + hl, Position.y + hl, Position.z - hl);


			//right
			glVertex3f(Position.x + hl, Position.y + hl, Position.z - hl);
			glVertex3f(Position.x + hl, Position.y - hl, Position.z - hl);
			glVertex3f(Position.x + hl, Position.y - hl, Position.z + hl);
			glVertex3f(Position.x + hl, Position.y + hl, Position.z + hl);


			//left
			glVertex3f(Position.x - hl, Position.y + hl, Position.z - hl);
			glVertex3f(Position.x - hl, Position.y - hl, Position.z - hl);
			glVertex3f(Position.x - hl, Position.y - hl, Position.z + hl);
			glVertex3f(Position.x - hl, Position.y + hl, Position.z + hl);

			//front

			glVertex3f(Position.x + hl, Position.y - hl, Position.z + hl);
			glVertex3f(Position.x - hl, Position.y - hl, Position.z + hl);
			glVertex3f(Position.x - hl, Position.y + hl, Position.z + hl);
			glVertex3f(Position.x + hl, Position.y + hl, Position.z + hl);

			//top

			glVertex3f(Position.x - hl, Position.y + hl, Position.z - hl);
			glVertex3f(Position.x - hl, Position.y + hl, Position.z + hl);
			glVertex3f(Position.x + hl, Position.y + hl, Position.z + hl);
			glVertex3f(Position.x + hl, Position.y + hl, Position.z - hl);

			glEnd();
		}
	}


};
rigid_cube* fire_cube;


class Node {
public:
	vec3 Position;
	vec3 Velocity = vec3(0.0f, 0.0f, 0.0f);
	vec3 Normal;
	vec3 force = vec3(0.0f, 0.0f, 0.0f);

	//for normal vector calculation
	std::vector <vec3> norms;

	bool isFixed;
	float mass = 1.0f;

	
	void add_force(vec3 additional_force) {
		force += additional_force;
	}

	//dt = timestamp calculate next position
	void integrate(double dt) {
		if (!isFixed) {
			vec3 a = force / mass;
			Velocity += a * (float)dt;
			
			Position += Velocity * (float)dt;
		}
		force.x = force.y = force.z = 0.0f;
	}

	Node(vec3 pos) {
		Position = pos;
	}

};
class mass_spring {
public:
	Node *N1;
	Node *N2;

	float spring_coef;
	float damping_coef = 0.05f;
	float naturalLength;

	//dt = timestep  ==> node의 force에 spring force 저장
		void internal_force(double dt) {

		vec3 direction = N1->Position - N2->Position;
		vec3 revel = N1->Velocity - N2->Velocity;
		float sLength = length(direction);

		//calculate hook's law, damping
		vec3 hook_force = (-1.0f)*spring_coef *(sLength - naturalLength)* direction / sLength;
		vec3 damping_force = (-1.0f)*damping_coef*revel;
		
		//add force
		N1->add_force(hook_force+damping_force);
		N2->add_force(-hook_force+damping_force);
	}

	mass_spring(Node* p1,  Node* p2) {
		
		N1 = p1;
		N2 = p2;
	}


};

class mass_cloth {
public:
	
	std::vector<Node *>nodes;
	std::vector<mass_spring *> spring;
	std::vector<Node *> faces;

	int size_x, size_y, size_z;
	double dx, dy, dz;
	double structural_coef;
	double shear_coef;
	double bending_coef;
	int iteration_n;


	vec3 startPos;

	int getX(int indx) {
		return (indx % size_x);
	}
	int getY(int indx) {
		return (indx / size_y);
	}


	void init(){
		nodes.clear();
		spring.clear();
		faces.clear();
		//node creation
		for (int i = 0; i < size_x; i++) {
			for (int j = 0; j < size_y; j++) {
				Node* xp = new Node(vec3(startPos.x + j, startPos.y , startPos.z+i));
				//fixed for four corners
				if ((j == 0 && i == 0) || (j == size_x - 1 && i == 0) || (j==0 && i ==size_y-1)||(j==size_x-1 &&i==size_y-1))
					xp->isFixed = true;
				else
					xp->isFixed = false;
				nodes.push_back(xp);
			}
		}

		//make structural spring Horizontal
		for (int i = 0; i < nodes.size(); i++) {
			if (getX(i) < size_x - 1) {
				mass_spring *sp1 = new mass_spring(nodes[i], nodes[i + 1]);
				sp1->spring_coef = structural_coef;
				sp1->naturalLength = dx;
				spring.push_back(sp1);
			}
		}

		//make structural spring Vertical
		for (int i = 0; i < nodes.size() - size_x; i++) {
			mass_spring *sp2 = new mass_spring(nodes[i], nodes[i + size_x]);
			sp2->spring_coef = structural_coef;
			sp2->naturalLength = dy;
			spring.push_back(sp2);
		}
		//make shear spring
		for (int i = 0; i < nodes.size() - size_x; i++) {
			if (getX(i) < size_x - 1) {
				mass_spring *sp3 = new mass_spring(nodes[i], nodes[i + size_x + 1]);
				sp3->spring_coef = shear_coef;
				sp3->naturalLength = sqrt(pow(dx, 2) + pow(dy, 2));
				spring.push_back(sp3);
			}
			if (getX(i) > 0) {
				mass_spring *sp4 = new mass_spring(nodes[i], nodes[i + size_x - 1]);
				sp4->spring_coef = shear_coef;
				sp4->naturalLength = sqrt(pow(dx, 2) + pow(dy, 2));
				spring.push_back(sp4);
			}

		}
		//make horizontal bending spring
		for (int i = 0; i < nodes.size(); i++) {
			if (getX(i) < size_x - 2) {
				mass_spring *sp5 = new mass_spring(nodes[i], nodes[i + 2]);
				sp5->spring_coef = bending_coef;
				sp5->naturalLength = 2*dx;
				spring.push_back(sp5);
			}
		}

		//make structural spring Vertical
		for (int i = 0; i < nodes.size() - (size_x*2); i++) {
			mass_spring *sp6 = new mass_spring(nodes[i], nodes[i + (2*size_x)]);
			sp6->spring_coef = bending_coef;
			sp6->naturalLength = 2*dy;
			spring.push_back(sp6);
		}

		//generate face
		for (int i = 0; i < nodes.size() - size_x; i++) {
			if (getX(i) < size_x - 1) {
				faces.push_back(nodes[i]);
				faces.push_back(nodes[i + 1]);
				faces.push_back(nodes[i + 1 + size_x]);

				faces.push_back(nodes[i]);
				faces.push_back(nodes[i + size_x]);
				faces.push_back(nodes[i + size_x + 1]);
			}
		}

	}

	void draw() {
		glEnable(GL_POLYGON_SMOOTH);
		glColor3f(1.0f, 0.0f, 0.0f);
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < faces.size(); i++) {
			//draw with faces and normal
			vec3 pos = faces[i]->Position;
			vec3 norm = faces[i]->Normal;
			glVertex3f(pos.x, pos.y, pos.z);
			glNormal3f(norm.x, norm.y, norm.z);
		}
		glEnd();
		//white line drawing
		if (strip) {
			glBegin(GL_LINE_STRIP);
			glColor3f(1.0f, 1.0f, 1.0f);
			for (int i = 0; i < faces.size(); i++) {
				vec3 pos = faces[i]->Position;
				glVertex3f(pos.x, pos.y, pos.z);

			}
			glEnd();
		}
	}
	// compute force(internal+external)
	void compute_force(double dt, vec3 gravity) {
		for (int i = 0; i < nodes.size(); i++) {
			//add gravity to all nodes
			nodes[i]->add_force(gravity*(nodes[i]->mass));
		}
		for (int i = 0; i < spring.size(); i++) {
			//add spring force to all springs
			spring[i]->internal_force(dt);
		}
	}

	void integrate(double dt) {
		//integrate forces -> get position
		for (int i = 0; i < nodes.size(); i++)
			nodes[i]->integrate(dt);
	}

	void computeNormal() {
		//normal calculation
		for (int i = 0; i < faces.size(); i += 3) {
			//compute face normal
			vec3 AB = faces[i + 1]->Position - faces[i]->Position;
			vec3 BC = faces[i + 2]->Position - faces[i + 1]->Position;
			vec3 Norm = normalize(cross(AB, BC));

			(faces[i]->norms).push_back(Norm);
			(faces[i + 1]->norms).push_back(Norm);
			(faces[i + 2]->norms).push_back(Norm);
		}

		//compute mean normal
		for (int i = 0; i < nodes.size(); i++) {
			vec3 sum(0.0f, 0.0f, 0.0f);
			for (int j = 0; j < (nodes[i]->norms).size(); j++) {
				sum += (nodes[i]->norms)[j];
			}
			nodes[i]->Normal = sum / (float)((nodes[i]->norms).size());
			nodes[i]->norms.clear();
		}
	}

	void collision() {
		int count = 0;
		vec3 min, max;
		//additional 1.5f -> corner bending
		float col_hl = fire_cube->hl +1.5;

		//range of top side of cube
		min = vec3((fire_cube->Position.x) - (col_hl), (fire_cube->Position.y) - (col_hl), (fire_cube->Position.z) - (col_hl));
		max = vec3((fire_cube->Position.x) + (col_hl), (fire_cube->Position.y) + (col_hl), (fire_cube->Position.z) + (col_hl));

		for (int i = 0; i < nodes.size(); i++) {
			vec3 pos = nodes[i]->Position;
			//collision detected between top of cube and cloth
			if ((pos.x >= min.x) && (pos.x <= max.x) && (pos.y >= min.y) && (pos.y <= max.y) && (pos.z >= min.z) && (pos.z <= max.z)) {
				//stop moving
				nodes[i]->Velocity *= 0.0f;
			}
			//collision detected between ground and cloth
			if (pos.y < ground.y) {
				//stop moving
				nodes[i]->Velocity *= 0.0f;
			}

			if (!isrender) { // count near nodes for firing
				if(length(pos - fire_cube->Position) <= (fire_cube->hl)*2)
					count++;
			}
		}
		if (count >= nodes.size()*0.8)
			isrender = true;
	}

	//unfix nodes
	void setFree() {
		for (int i = 0; i < nodes.size(); i++) {
			if (nodes[i]->isFixed)
				nodes[i]->isFixed = false;
		}
	}

};
mass_cloth* cloth;

class Simulator
{
public:

	void drawGround(void) {
		glEnable(GL_BLEND);
		glEnable(GL_TEXTURE_2D);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);

		glBindTexture(GL_TEXTURE_2D, mTexture[1]);

		glBegin(GL_QUADS);
		glColor3f(1.0f, 1.0f, 1.0f);
		glNormal3f(0.0, 1.0, 0.0);
		glTexCoord2f(0.0f, 0.0f);
		glVertex3f(-LENGTH, -10.0f-(LENGTH/2), -LENGTH);
		glTexCoord2f(0.0f, 1.0f);
		glVertex3f(LENGTH, -10.0f - (LENGTH / 2), -LENGTH);
		glTexCoord2f(1.0f, 1.0f);
		glVertex3f(LENGTH, -10.0f - (LENGTH / 2), LENGTH);
		glTexCoord2f(1.0f, 0.0f);
		glVertex3f(-LENGTH, -10.0f - (LENGTH / 2), LENGTH);
		glEnd();

		glDisable(GL_BLEND);
		glDisable(GL_TEXTURE_2D);

	}
	void Initialize(void)
	{
		timeStep = 0.01;
		ParticleSystem.init(NUMBER_OF_FIRE);

	}
	void C_Initialize(void) {

		cloth = new mass_cloth();
		cloth->dx = 1;
		cloth->dy = 1;
		cloth->dz = 1;
		cloth->size_x = 50;
		cloth->size_y = 50;
		cloth->structural_coef = 1000;
		cloth->shear_coef = 50;
		cloth->bending_coef = 1000;
		cloth->iteration_n = 20;

		cloth->startPos = vec3(-25.0f,(LENGTH/2), -25.0f);
		cloth->init();
	}

	//initialize cube
	void R_Initialize(void) {
		fire_cube = new rigid_cube(vec3(0,  -(LENGTH / 2),0), 10);
		fire_cube->Color = vec3(1.0f, 0.5f, 0.0f);
	}

	void Lighting() {
		glShadeModel(GL_SMOOTH);

		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_COLOR_MATERIAL);

		float light_ambient[] = { 0.3, 0.3, 0.3, 1.0 };
		float light_diffuse[] = { 0.8, 0.8, 0.8, 1.0 };
		float light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		float light_pos[] = { 0.0f, 150.0f, 0.0f, 0.0f };

		glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
		glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
		
	}
	void Update()
	{
		//iteration for iteration_n
		vec3 gravity(0.0, -9.8 / cloth->iteration_n, 0.0);
		for (int iter = 0; iter < cloth->iteration_n; iter++) {
			
			cloth->compute_force(timeStep, gravity);
			cloth->integrate(timeStep);
			//collision detection
			cloth->collision();
		}
		cloth->computeNormal();
		
		//available when isrender and master is true
		if (isrender && master) {
			ParticleSystem.Movement(timeStep);
			ParticleSystem.Coloring();
		}

	}

	void Render() {
		drawGround();
		fire_cube->draw();
		cloth->draw();

		//rendering  when isrender and master is true
		if (isrender &&master) {
			ParticleSystem.draw();
		}
	}

};
Simulator S_Simulator;


class Viewer
{
public:

	void Initialize(void) {
		LoadMeshFromFile("fire.bmp", 0);
		LoadMeshFromFile("floor.bmp", 1);
		S_Simulator.Initialize();
		S_Simulator.C_Initialize();
		S_Simulator.R_Initialize();
		S_Simulator.Lighting();
	}
	void Update() {

			S_Simulator.Update();
	}

	void Render(void) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(rotate[0], 1.0, 0.0, 0.0);
		glRotatef(rotate[1], 0.0, 1.0, 0.0);
		glRotatef(rotate[2], 0.0, 0.0, 1.0);
		glOrtho(-(LENGTH + 20), (LENGTH + 20), -(LENGTH + 20), (LENGTH + 20), -(LENGTH + 20), (LENGTH + 20));
		glPushMatrix();

		S_Simulator.Render();
		glPopMatrix();

		glutSwapBuffers();
		glutPostRedisplay();


	}

	void Keyboard(unsigned char key, int x, int y) {
		switch (key) {
		//falling cloth
		case ' ':
			cloth->setFree();
			break;

		//rotate
		case 'a':
		case 'A':
			rotate[0] -= 10;
			break;

		case 's':
		case 'S':
			rotate[1] -= 10;
			break;

		case 'd':
		case 'D':
			rotate[2] -= 10;
			break;

		case 'r': //reset
		case 'R':
			rotate[0] =rotate[1] =rotate[2]= 0;
			break;
		
		case 'f'://fire on/off
		case 'F':
			isrender = !isrender;
			break;
		
		case 'n'://reset simulation
		case 'N':
			isrender = false;
			master = true;
			Initialize();
			break;

		
		case 'l'://draw white line
		case 'L':
			strip = !strip;
			break;

		case '0': //hide fire
			master = !master;
			break;
		}
	}

	void Reshape(int w, int h) {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(45, (float)2 / h, 0.1, 500);
		glViewport(0, 0, w, h);
	}

};
Viewer OpenGL_Viewer;



void Initialize(void) {
	OpenGL_Viewer.Initialize();
}
void Update(int value) {
	OpenGL_Viewer.Update();
	glutPostRedisplay();
	glutTimerFunc(10, Update, 0);
}
void Render(void) {
	OpenGL_Viewer.Render();
}
void Keyboard(unsigned char key, int x, int y) {
	OpenGL_Viewer.Keyboard(key, x, y);
}

void Reshape(int w, int h) {
	OpenGL_Viewer.Reshape(w, h);
}

int main(int argc, char* argv[]) {
	glutInitWindowSize(800, 700);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);

	glEnable(GL_DEPTH_TEST);

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	glutCreateWindow("TermProject");
	rotate[0] = -20;
	rotate[1] = -50;

	Initialize();
	glutTimerFunc(10, Update, 0);
	glutDisplayFunc(Render);
	
	glutKeyboardFunc(Keyboard);

	//glutReshapeFunc(Reshape);

	glutMainLoop();

	return 0;
}