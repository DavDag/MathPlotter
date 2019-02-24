#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include "exp/mathexpr.h"
#include <bits/stdc++.h>

#define CENTERED 0
#define ALIGN_RIGHT 1

const double pi = std::atan(1.0f) * 4;
const double e = std::exp(1.0f);
const double inf = std::numeric_limits<double>::infinity();

GLFWwindow *window;
unsigned int vertexbuffer, axisbuffer, numbersbuffer, gridbuffer, deepgridbuffer, integratebuffer, guidelinesbuffer, parabolasbuffer;

const int num_points = 1000, intervalError = 50000;
double scale = 0.1f, xoffset = 0, yoffset = 0, stepstart = 0, stepend = 0, area = 0.0f, errore = 0.0f;
int stepintegral = 1, stepforcheckinterval = 1000, parabolaprecision = 100, intervalMinMaxCheckPrecision = 1000;
double points[num_points << 1];
std::vector<double> axis, numbers, grid, deepgrid, integrals, guidelines, parabolas;
bool precision_scale = true, drag = false, integrable = false, check = false, updatefunction = true;
int integrate = -1;
bool breakable = false;
double screenWidth = 1100.0f, screenHeight = 600.0f;
std::pair<double, double> lastMousePos = {0.0f, 0.0f}, currMousePos = {0.0f, 0.0f};

char function[500] = "x";
double funx;
RVar *xvar;
RVar* vararray[1];
ROperation *op;

void MousePosProc(GLFWwindow*, double, double);
void MouseScrollProc(GLFWwindow*, double, double);
void MouseButtonProc(GLFWwindow*, int, int, int);
void KeyboardProc(GLFWwindow*, int, int, int, int);
void APIENTRY ErrorProc(GLenum, GLenum, GLuint, GLenum, GLsizei, const GLchar*, const void*);
void Init(int, int, const std::string&);
void Destroy();
void Update();
void DrawContent();
void DrawGui();
void ElaborateFunction();
void CheckForInterval();
void CalcError();
void DrawParabola(double, double, double);
void DrawIntegral();
void DrawGrid(double, double, std::string);
std::vector<double> DrawNumber(double, double, unsigned int, double);
double Func(double);
double FirstDeriv(double);
double SecondDeriv(double);
double ThirdDeriv(double);
double FourthDeriv(double);

int main()
{
	int scelta = 0;
	do
	{
		std::cout << "Insert Screen Size : ";
		std::cin >> screenWidth >> screenHeight;
		Init(screenWidth, screenHeight, "Grafi");
		while (true && !glfwWindowShouldClose(window))
	    {
	    	Update();
	    	DrawContent();
			DrawGui();
			glfwSwapBuffers(window);
			glfwPollEvents();
	    }
	    Destroy();
	    std::cout << "Premere 1 per continuare : ";
	    std::cin >> scelta;
    } while(scelta == 1);
	return EXIT_SUCCESS;
}

double FirstDeriv(double x)
{
	double h = 0.0001f;
	double x2 = x + h;
	double x1 = x - h;
	double f2 = Func(x2);
	double f1 = Func(x1);
	double approxDerivative = (f2 - f1) / (2.0 * h);
	return approxDerivative;
}

double SecondDeriv(double x)
{
	double h = 0.0001f;
	double x2 = x + h;
	double x1 = x - h;
	double f2 = FirstDeriv(x2);
	double f1 = FirstDeriv(x1);
	double approxDerivative = (f2 - f1) / (2.0 * h);
	return approxDerivative;
}

double ThirdDeriv(double x)
{
	double h = 0.0001f;
	double x2 = x + h;
	double x1 = x - h;
	double f2 = SecondDeriv(x2);
	double f1 = SecondDeriv(x1);
	double approxDerivative = (f2 - f1) / (2.0 * h);
	return approxDerivative;
}

double FourthDeriv(double x)
{
	double h = 0.0001f;
	double x2 = x + h;
	double x1 = x - h;
	double f2 = ThirdDeriv(x2);
	double f1 = ThirdDeriv(x1);
	double approxDerivative = (f2 - f1) / (2.0 * h);
	return approxDerivative;
}

double Func(double x)
{
	double y = 0;
	funx = x;
	y = op->Val();
	if(y == inf && breakable) integrable = false;
	return y;
}

void ElaborateFunction()
{
	if(!updatefunction) return;
	updatefunction = false;
	xvar = new RVar("x", &funx);
	vararray[0] = xvar;
	op = new ROperation(function, 1, vararray);
}

void Update()
{
	// Disegno la funzione
	ElaborateFunction();
	double x = (-screenWidth / 2 + xoffset / scale) * scale;
	double func_step = (screenWidth * scale) / ((double)num_points + 1);
	for(int i = 0; i < num_points; ++i)
	{
		points[i << 1] = (x - xoffset) / scale / (screenWidth / 2);
		points[(i << 1) + 1] = (Func(x) - yoffset) / scale / (screenHeight / 2);
		x += func_step;
	}
	// Disegno gli assi cartesiani
	axis = {
		 1.0f, -yoffset / scale / screenHeight * 2.0f,
		-1.0f, -yoffset / scale / screenHeight * 2.0f,
		-xoffset / scale / screenWidth * 2.0f,  1.0f,
		-xoffset / scale / screenWidth * 2.0f, -1.0f
	};
	// Disegno i numeri lungo gli assi
	numbers.clear();
	grid.clear();
	deepgrid.clear();
	double step = scale * 100;
	// Ascisse positive
	for(double x = -xoffset; x < screenWidth / 2 * scale; x += step)
	{
		std::vector<double> num = DrawNumber(x + xoffset, 8.75f, CENTERED, scale);
		for(int i = 0; i < num.size(); i += 2) num[i] = (num[i] + x / scale) / (screenWidth / 2);
		for(int i = 1; i < num.size(); i += 2) num[i] = (num[i] - yoffset / scale) / (screenHeight / 2);
		for(auto el : num) numbers.push_back(el);
		// Disegno la griglia
		DrawGrid(x, x + step, "x");
	}
	// Ascisse negative
	DrawGrid(-xoffset, -xoffset - step, "x");
	for(double x = -xoffset - step; x > - screenWidth / 2 * scale; x -= step)
	{
		std::vector<double> num = DrawNumber(x + xoffset, 8.75f, CENTERED, scale);
		for(int i = 0; i < num.size(); i += 2) num[i] = (num[i] + x / scale) / (screenWidth / 2);
		for(int i = 1; i < num.size(); i += 2) num[i] = (num[i] - yoffset / scale) / (screenHeight / 2);
		for(auto el : num) numbers.push_back(el);
		// Disegno la griglia
		DrawGrid(x, x - step, "x");
	}
	// Ordinate positive
	DrawGrid(-yoffset, -yoffset + step, "y");
	for(double y = -yoffset + step; y < screenHeight / 2 * scale; y += step)
	{
		std::vector<double> num = DrawNumber(y + yoffset, 8.75f, ALIGN_RIGHT, scale);
		for(int i = 0; i < num.size(); i += 2) num[i] = (num[i] - xoffset / scale) / (screenWidth / 2);
		for(int i = 1; i < num.size(); i += 2) num[i] = (num[i] + y / scale) / (screenHeight / 2);
		for(auto el : num) numbers.push_back(el);
		// Disegno la griglia
		DrawGrid(y, y + step, "y");
	}
	// Ordinate negative
	DrawGrid(-yoffset, -yoffset - step, "y");
	for(double y = -yoffset - step; y > -screenHeight / 2 * scale; y -= step)
	{
		std::vector<double> num = DrawNumber(y + yoffset, 8.75f, ALIGN_RIGHT, scale);
		for(int i = 0; i < num.size(); i += 2) num[i] = (num[i] - xoffset / scale) / (screenWidth / 2);
		for(int i = 1; i < num.size(); i += 2) num[i] = (num[i] + y / scale) / (screenHeight / 2);
		for(auto el : num) numbers.push_back(el);
		// Disegno la griglia
		DrawGrid(y, y - step, "y");
	}
	// Disegno x e F(x) con x = posizione del mouse
	x = (currMousePos.first - screenWidth / 2) * scale;
	double y = Func(x);
	axis.push_back((x) / scale / (screenWidth / 2));
	axis.push_back((-yoffset) / scale / (screenHeight / 2));
	axis.push_back((x) / scale / (screenWidth / 2));
	axis.push_back((Func(x + xoffset) - yoffset) / scale / (screenHeight / 2));
	axis.push_back((-xoffset) / scale / (screenWidth / 2));
	axis.push_back((Func(x + xoffset) - yoffset) / scale / (screenHeight / 2));
	axis.push_back((x) / scale / (screenWidth / 2));
	axis.push_back((Func(x + xoffset) - yoffset) / scale / (screenHeight / 2));
	// Disegno i numeri sopra il mouse
	/*x += xoffset;
	y -= yoffset;
	double xshift = 9 * std::log10(std::abs(x)) + 18, yshift = (y < 0) ? -10 : 10;
	x = (double)((int)(x * 100)) / 100.0f;
	y = (double)((int)(y * 100)) / 100.0f;
	std::vector<double> num = drawNumber(x, 8.75f, CENTERED, scale);
	for(int i = 0; i < num.size(); i += 2) num[i] = (num[i] + x / scale - xshift) / (screenWidth / 2);
	for(int i = 1; i < num.size(); i += 2) num[i] = (num[i] + y / scale + yshift) / (screenHeight / 2);
	for(auto el : num) numbers.push_back(el);
	num = drawNumber(y, 8.75f, CENTERED, scale);
	for(int i = 0; i < num.size(); i += 2) num[i] = (num[i] + x / scale + xshift) / (screenWidth / 2);
	for(int i = 1; i < num.size(); i += 2) num[i] = (num[i] + y / scale + yshift) / (screenHeight / 2);
	for(auto el : num) numbers.push_back(el);*/
	// Disegno l'integrale
	breakable = true;
	DrawIntegral();
	breakable = false;
}

void CalcError()
{
	errore = 0.0f;
	double step = (stepend - stepstart) / (double) intervalError;
	double M = 0.0f;
	for(double i = stepstart; i < stepend; i += step)
	{
		if(integrate == 0 || integrate == 1) M = std::max(M, std::abs(FirstDeriv(i)));
		else if(integrate == 2) M = std::max(M, std::abs(SecondDeriv(i)));
		else if(integrate == 3) M = std::max(M, std::abs(FourthDeriv(i)));
	}
	if(integrate == 0 || integrate == 1) errore = std::pow((stepend - stepstart), 2) / (2.0f * (double) stepintegral) * M;
	else if(integrate == 2) errore = std::pow((stepend - stepstart), 3) / (12.0f * std::pow((double) stepintegral, 2)) * M;
	else if(integrate == 3) errore = std::pow((stepend - stepstart), 5) / (2880.0f * std::pow((double) stepintegral, 4)) * M;
}

void CheckForInterval()
{
	if(!check) return;
	check = false;
	double step = (stepend - stepstart) / (double)stepforcheckinterval;
	int sign = 0;
	bool allgood = true;
	for(double i = stepstart; i < stepend; i += step)
	{
		double y = Func(i);
		if(y == 0 || y == -0) continue;
		if(sign == 0) sign = (y > 0) ? 1 : -1;
		else if(sign == -1 && y < 0) continue;
		else if(sign == +1 && y > 0) continue;
		else
		{
			allgood = false;
			break;
		}
	}
	integrable = allgood;
}

void DrawParabola(double x1, double x2, double x3)
{
	double y1 = Func(x1), y2 = Func(x2), y3 = Func(x3);
//	printf("[%f, %f], [%f, %f], [%f, %f]\n", x1, y1, x2, y2, x3, y3);
	// Per calcolare i coefficenti uso i "polinomi di Legendre"
	// credit http://www.francococca.com/gg/parabola_per_tre_punti.asp
	double a = (y1)/((x1-x2)*(x1-x3))+(y2)/((x2-x1)*(x2-x3))+(y3)/((x3-x1)*(x3-x2));
	double b = -(y1*(x2+x3)/((x1-x2)*(x1-x3))+y2*(x1+x3)/((x2-x1)*(x2-x3))+y3*(x1+x2)/((x3-x1)*(x3-x2)));
	double c = y1*x2*(x3)/((x1-x2)*(x1-x3))+y2*x1*(x3)/((x2-x1)*(x2-x3))+y3*x1*(x2)/((x3-x1)*(x3-x2));
	std::string str_func = std::to_string(a) + "x^2 + " + std::to_string(b) + "x + " + std::to_string(c);
	char temp_func[str_func.size() + 1];
	std::strcpy(temp_func, str_func.c_str());
	double temp_funx;
	RVar temp_xvar("x", &temp_funx);
	RVar* temp_vararray[1];
	temp_vararray[0] = &temp_xvar;
	ROperation temp_op(temp_func, 1, temp_vararray);
//	std::cout << temp_op.Expr() << "\n";
	double x = x1;
	double func_step = ((x3 - x1)) / ((double)parabolaprecision);
	for(int i = 0; i <= parabolaprecision; ++i)
	{
		temp_funx = x;
		parabolas.push_back((temp_funx - xoffset) / scale / (screenWidth / 2));
		parabolas.push_back((temp_op.Val() - yoffset) / scale / (screenHeight / 2));
		x += func_step;
	}
	double delta = (x3 - x1) / (double) parabolaprecision;
	for(int i = 0; i <= parabolaprecision; ++i)
	{
		temp_funx = (x1 + delta * (i + 1));
		double massimo = temp_op.Val();
		temp_funx = (x1 + delta * i);
		double minimo = temp_op.Val();
		if(massimo < minimo) std::swap(minimo, massimo);

		integrals.push_back((x1 + delta * i - xoffset) / scale / (screenWidth / 2));
		integrals.push_back((-yoffset) / scale / (screenHeight / 2));
		integrals.push_back((x1 + delta * i - xoffset) / scale / (screenWidth / 2));
		temp_funx = (x1 + delta * i);
		integrals.push_back((temp_op.Val() - yoffset) / scale / (screenHeight / 2));
		integrals.push_back((x1 + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
		temp_funx = (x1 + delta * (i + 1));
		integrals.push_back((temp_op.Val() - yoffset) / scale / (screenHeight / 2));
	
//		guidelines.push_back((x1 + delta * i - xoffset) / scale / (screenWidth / 2));
//		guidelines.push_back((-yoffset) / scale / (screenHeight / 2));
//		guidelines.push_back((x1 + delta * i - xoffset) / scale / (screenWidth / 2));
//		temp_funx = (x1 + delta * i);
//		guidelines.push_back((temp_op.Val() - yoffset) / scale / (screenHeight / 2));
		
		integrals.push_back((x1 + delta * i - xoffset) / scale / (screenWidth / 2));
		integrals.push_back((-yoffset) / scale / (screenHeight / 2));
		integrals.push_back((x1 + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
		temp_funx = (x1 + delta * (i + 1));
		integrals.push_back((temp_op.Val() - yoffset) / scale / (screenHeight / 2));
		integrals.push_back((x1 + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
		integrals.push_back((-yoffset) / scale / (screenHeight / 2));
		
//		guidelines.push_back((x1 + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
//		guidelines.push_back((-yoffset) / scale / (screenHeight / 2));
//		guidelines.push_back((x1 + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
//		temp_funx = (x1 + delta * (i + 1));
//		guidelines.push_back((temp_op.Val() - yoffset) / scale / (screenHeight / 2));
	}
	guidelines.push_back((x1 - xoffset) / scale / (screenWidth / 2));
	guidelines.push_back((-yoffset) / scale / (screenHeight / 2));
	guidelines.push_back((x1 - xoffset) / scale / (screenWidth / 2));
	temp_funx = x1;
	guidelines.push_back((temp_op.Val() - yoffset) / scale / (screenHeight / 2));
	guidelines.push_back((x3 - xoffset) / scale / (screenWidth / 2));
	guidelines.push_back((-yoffset) / scale / (screenHeight / 2));
	guidelines.push_back((x3 - xoffset) / scale / (screenWidth / 2));
	temp_funx = x3;
	guidelines.push_back((temp_op.Val() - yoffset) / scale / (screenHeight / 2));
}

void DrawIntegral()
{
	integrals.clear();
	guidelines.clear();
	parabolas.clear();
	area = 0.0f;
	if(stepstart == stepend || integrate == -1) return;
	CheckForInterval();
	if(!integrable) return;
	if(integrate == 0 || integrate == 1 || integrate == 2)
	{
		double delta = (stepend - stepstart) / (double) stepintegral;
		for(int i = 0; i < stepintegral; ++i)
		{
			double massimo = Func(stepstart + delta * (i + 1)), minimo = Func(stepstart + delta * i);
			double tempx = stepstart + delta * (i);
			double tempstep = (delta) / (double) intervalMinMaxCheckPrecision;
			if(massimo < minimo) std::swap(massimo, minimo);
			for(int j = 0; j < intervalMinMaxCheckPrecision && integrate != 2; ++j)
			{
				double tempy = Func(tempx + j * tempstep);
				if(tempy < minimo) minimo = tempy;
				if(tempy > massimo) massimo = tempy;
			}
			
			if(integrate == 0) area += minimo * delta;
			else if(integrate == 1) area += massimo * delta;
			else area += (minimo + massimo) / 2.0f * delta;
			
			integrals.push_back((stepstart + delta * i - xoffset) / scale / (screenWidth / 2));
			integrals.push_back((-yoffset) / scale / (screenHeight / 2));
			integrals.push_back((stepstart + delta * i - xoffset) / scale / (screenWidth / 2));
			if(integrate == 0)		integrals.push_back((minimo - yoffset) / scale / (screenHeight / 2));
			else if(integrate == 1)	integrals.push_back((massimo - yoffset) / scale / (screenHeight / 2));
			else 					integrals.push_back((Func(stepstart + delta * i) - yoffset) / scale / (screenHeight / 2));
			integrals.push_back((stepstart + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
			if(integrate == 0)		integrals.push_back((minimo - yoffset) / scale / (screenHeight / 2));
			else if(integrate == 1)	integrals.push_back((massimo - yoffset) / scale / (screenHeight / 2));
			else 					integrals.push_back((Func(stepstart + delta * (i + 1)) - yoffset) / scale / (screenHeight / 2));
		
			guidelines.push_back((stepstart + delta * i - xoffset) / scale / (screenWidth / 2));
			guidelines.push_back((-yoffset) / scale / (screenHeight / 2));
			guidelines.push_back((stepstart + delta * i - xoffset) / scale / (screenWidth / 2));
			if(integrate == 0)		guidelines.push_back((minimo - yoffset) / scale / (screenHeight / 2));
			else if(integrate == 1) guidelines.push_back((massimo - yoffset) / scale / (screenHeight / 2));
			else 					guidelines.push_back((Func(stepstart + delta * i) - yoffset) / scale / (screenHeight / 2));
			
			integrals.push_back((stepstart + delta * i - xoffset) / scale / (screenWidth / 2));
			integrals.push_back((-yoffset) / scale / (screenHeight / 2));
			integrals.push_back((stepstart + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
			if(integrate == 0)		integrals.push_back((minimo - yoffset) / scale / (screenHeight / 2));
			else if(integrate == 1) integrals.push_back((massimo - yoffset) / scale / (screenHeight / 2));
			else 					integrals.push_back((Func(stepstart + delta * (i + 1)) - yoffset) / scale / (screenHeight / 2));
			integrals.push_back((stepstart + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
			integrals.push_back((-yoffset) / scale / (screenHeight / 2));
			
			guidelines.push_back((stepstart + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
			guidelines.push_back((-yoffset) / scale / (screenHeight / 2));
			guidelines.push_back((stepstart + delta * (i + 1) - xoffset) / scale / (screenWidth / 2));
			if(integrate == 0)		guidelines.push_back((minimo - yoffset) / scale / (screenHeight / 2));
			else if(integrate == 1)	guidelines.push_back((massimo - yoffset) / scale / (screenHeight / 2));
			else 					guidelines.push_back((Func(stepstart + delta * (i + 1)) - yoffset) / scale / (screenHeight / 2));
		}
	}
	else if(integrate == 3)
	{
		double delta = (stepend - stepstart) / ((double) stepintegral);
		for(int i = 0; i < stepintegral; ++i)
		{
			double x1 = stepstart + delta * i;
			double x3 = stepstart + delta * (i + 1);
			DrawParabola(x1, (x1 + x3) / 2.0f, x3);
			area += (Func(x1) + Func(x3) + 4 * Func((x1 + x3) / 2.0f)) * (delta / 2.0f / 3.0f);
		}
	}
}

void DrawGrid(double val, double dval, std::string axe)
{
	if(axe == "x")
	{
		double x = val;
		grid.push_back(x / scale / screenWidth * 2);
		grid.push_back(-1.0f);
		grid.push_back(x / scale / screenWidth * 2);
		grid.push_back(1.0f);
		for(int i = 1; i < 5; ++i)
		{
			x = (dval - val ) * i / 5.0f + val;
			deepgrid.push_back(x / scale / screenWidth * 2);
			deepgrid.push_back(-1.0f);
			deepgrid.push_back(x / scale / screenWidth * 2);
			deepgrid.push_back(1.0f);
		}
	}
	else if(axe == "y")
	{
		double y = val;
		grid.push_back(-1.0f);
		grid.push_back(y / scale / screenHeight * 2);
		grid.push_back(1.0f);
		grid.push_back(y / scale / screenHeight * 2);
		for(int i = 1; i < 5; ++i)
		{
			y = (dval - val ) * i / 5.0f + val;
			deepgrid.push_back(-1.0f);
			deepgrid.push_back(y / scale / screenHeight * 2);
			deepgrid.push_back(1.0f);
			deepgrid.push_back(y / scale / screenHeight * 2);
		}
	}
}

void DrawContent()
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	
	// Disegno la griglia ridotta
	glLineWidth(0.5f);
	glColor3f(0.9f, 0.9f, 0.9f);
	glBindBuffer(GL_ARRAY_BUFFER, deepgridbuffer);
	glBufferData(GL_ARRAY_BUFFER, deepgrid.size() * sizeof(double), deepgrid.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), 0);
	glDrawArrays(GL_LINES, 0, deepgrid.size() >> 1);
	
	// Disegno la griglia
	glLineWidth(0.5f);
	glColor3f(0.75f, 0.75f, 0.75f);
	glBindBuffer(GL_ARRAY_BUFFER, gridbuffer);
	glBufferData(GL_ARRAY_BUFFER, grid.size() * sizeof(double), grid.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), 0);
	glDrawArrays(GL_LINES, 0, grid.size() >> 1);
	
	// Disegno l'integrazione
	glColor3f(93.0f / 255.0f, 173.0f / 255.0f, 226.0f / 255.0f);
	glBindBuffer(GL_ARRAY_BUFFER, integratebuffer);
	glBufferData(GL_ARRAY_BUFFER, integrals.size() * sizeof(double), integrals.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), 0);
	glDrawArrays(GL_TRIANGLES, 0, integrals.size() >> 1);
	
	// Disegno la funzione
	glLineWidth(1.0f);
	glColor3f(0.2f, 0.0f, 1.0f);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(points), points, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), 0);
	glDrawArrays(GL_LINE_STRIP, 0, num_points);
	
	// Disegno le parabole
	glLineWidth(1.5f);
	glColor3f(1.0f, 0.0f, 0.0f);
	glBindBuffer(GL_ARRAY_BUFFER, parabolasbuffer);
	glBufferData(GL_ARRAY_BUFFER, parabolas.size() * sizeof(double), parabolas.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), 0);
	glDrawArrays(GL_LINE_STRIP, 0, parabolas.size() >> 1);
	
	// Disegno le linee guida
	glLineWidth(1.0f);
	glColor3f(91.0f / 255.0f, 84.0f / 255.0f, 111.0f / 255.0f);
	glBindBuffer(GL_ARRAY_BUFFER, guidelinesbuffer);
	glBufferData(GL_ARRAY_BUFFER, guidelines.size() * sizeof(double), guidelines.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), 0);
	glDrawArrays(GL_LINES, 0, guidelines.size() >> 1);
	
	// Disegno gli assi
	glLineWidth(1.0f);
	glColor3f(0.3f, 0.3f, 0.3f);
	glBindBuffer(GL_ARRAY_BUFFER, axisbuffer);
	glBufferData(GL_ARRAY_BUFFER, axis.size() * sizeof(double), axis.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), 0);
	glDrawArrays(GL_LINES, 0, axis.size() >> 1);
		
	// Disegno i valori numerici
	glLineWidth(0.5f);
	glColor3f(0.0f, 0.0f, 0.0f);
	glBindBuffer(GL_ARRAY_BUFFER, numbersbuffer);
	glBufferData(GL_ARRAY_BUFFER, numbers.size() * sizeof(double), numbers.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), 0);
	glDrawArrays(GL_LINES, 0, numbers.size() >> 1);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void DrawGui()
{
	ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    ImGui::Begin("Integrazione Numerica");
	ImGui::Text("F(x): "); ImGui::SameLine();
	if(ImGui::InputText("", function, 500)) updatefunction = check = true;
	//ImGui::Text("Scala: %.3f", scale);
	if(!integrable) ImGui::Text("Inserire un intervallo in cui la\nfunzione e': \n -STRETTAMENTE CRESCENTE O DECRESCENTE \n -CONTINUA");
	if(ImGui::InputDouble(": a", &stepstart)) check = true;
	if(ImGui::InputDouble(": b", &stepend)) check = true;
	if(ImGui::InputInt(": n", &stepintegral)) if(stepintegral > 1) CalcError();
	if(integrable)
	{
		ImGui::Text("Area = %lf", ((area >= 0) ? area : -area));
		ImGui::Text("Errore = %lf", errore);
	}
	if(ImGui::RadioButton("Nessuna", &integrate, -1)) CalcError();
	if(ImGui::RadioButton("Rettangoli per Difetto", &integrate, 0)) CalcError();
	if(ImGui::RadioButton("Rettangoli per Eccesso", &integrate, 1)) CalcError();
	if(ImGui::RadioButton("Trapezi", &integrate, 2)) CalcError();
	if(ImGui::RadioButton("Parabole", &integrate, 3)) CalcError();
	ImGui::End();
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
	if(scale < 0.001f) scale = 0.001f;
	else if(scale > 1000.0f) scale = 1000.0f;
	if(stepend < stepstart) stepend = stepstart;
	if(stepintegral < 1) stepintegral = 1;
	if(stepintegral > 500) stepintegral = 500;
}

std::vector<double> nums[] = {
	{ // Zero
		 0.2,  0.2,	 0.2,  0.8,
		 0.2,  0.8,  0.3, 0.95,
		 0.3, 0.95,  0.4,  1.0,
		 0.4,  1.0,  0.6,  1.0,
		 0.6,  1.0,  0.7, 0.95,
		 0.7, 0.95,  0.8,  0.8,
		 0.8,  0.8,  0.8,  0.2,
		 0.8,  0.2,  0.7, 0.05,
		 0.7, 0.05,  0.6,  0.0,
		 0.6,  0.0,  0.4,  0.0,
		 0.4,  0.0,  0.3, 0.05,
		 0.3, 0.05,  0.2,  0.2,
	},
	{ // Uno
		 0.2,  0.7,  0.8,  1.0,
		 0.8,  1.0,  0.8,    0
	},
	{ // Due
		 0.3,  0.6,  0.2,  0.8,
		 0.2,  0.8,  0.4,  1.0,
		 0.4,  1.0,  0.6,  1.0,
		 0.6,  1.0,  0.8,  0.8,
		 0.8,  0.8,  0.8,  0.6,
		 0.8,  0.6,  0.2,  0.0,
		 0.2,  0.0,  0.8,  0.0
	},
	{ // Tre
		 0.2,  1.0,  0.8,  1.0,
		 0.8,  1.0,  0.5,  0.6,
		 0.5,  0.6,  0.7,  0.6,
		 0.7,  0.6,  0.8,  0.5,
		 0.8,  0.5,  0.8,  0.2,
		 0.8,  0.2,  0.6,  0.0,
		 0.6,  0.0,  0.3,  0.0,
		 0.3,  0.0,  0.2,  0.1
	},
	{ // Quattro
		 0.6,  0.0,  0.6,  1.0,
		 0.6,  1.0,  0.2,  0.3,
		 0.2,  0.3,  0.8,  0.3
	},
	{ // Cinque
		 0.2,  0.1,  0.4,  0.0,
		 0.4,  0.0,  0.6,  0.0,
		 0.6,  0.0,  0.8,  0.2,
		 0.8,  0.2,  0.8,  0.5,
		 0.8,  0.5,  0.6,  0.7,
		 0.6,  0.7,  0.4,  0.7,
		 0.4,  0.7,  0.2,  0.6,
		 0.2,  0.6,  0.2,  1.0,
		 0.2,  1.0,  0.8,  1.0
	},
	{ // Sei
		 0.8,  1.0,  0.4,  0.9,
		 0.4,  0.9,  0.2,  0.6,
		 0.2,  0.6,  0.2,  0.2,
		 0.2,  0.2,  0.4,  0.0,
		 0.4,  0.0,  0.6,  0.0,
		 0.6,  0.0,  0.8,  0.2,
		 0.8,  0.2,  0.8,  0.4,
		 0.8,  0.4,  0.6,  0.5,
		 0.6,  0.5,  0.4,  0.5,
		 0.4,  0.5,  0.2,  0.4
	},
	{ // Sette
		 0.2,  1.0,  0.8,  1.0,
		 0.8,  1.0,  0.2,  0.0
	},
	{ // Otto
		 0.4,  0.0,  0.2, 0.15,
		 0.2, 0.15,  0.2,  0.4,
		 0.2,  0.4,  0.4,  0.5,
		 0.4,  0.5,  0.2,  0.6,
		 0.2,  0.6,  0.2, 0.85,
		 0.2, 0.85,  0.4,  1.0,
		 0.4,  1.0,  0.6,  1.0,
		 0.6,  1.0,  0.8, 0.85,
		 0.8, 0.85,  0.8,  0.6,
		 0.8,  0.6,  0.6,  0.5,
		 0.6,  0.5,  0.8,  0.4,
		 0.8,  0.4,  0.8, 0.15,
		 0.8, 0.15,  0.6,  0.0,
		 0.6,  0.0,  0.4,  0.0,
		 
		 0.4,  0.5,  0.6,  0.5
	},
	{ // Nove
		 0.2,  0.0,  0.6,  0.1,
		 0.6,  0.1,  0.8,  0.4,
		 0.8,  0.4,  0.8,  0.8,
		 0.8,  0.8,  0.6,  1.0,
		 0.6,  1.0,  0.4,  1.0,
		 0.4,  1.0,  0.2,  0.8,
		 0.2,  0.8,  0.2,  0.6,
		 0.2,  0.6,  0.4,  0.5,
		 0.4,  0.5,  0.6,  0.5,
		 0.6,  0.5,  0.8,  0.6
	},
	{ // Virgola
		 0.5,  0.0,  0.2, -0.4
	},
	{ // Meno
		 0.2,  0.5,  0.8,  0.5
	}
	
};

std::vector<double> DrawNumber(double number, double scale, unsigned int option, double worldscale)
{
	double PRECISION = 10.0f;
	std::vector<int> num(1, 0);
	if(worldscale < 1.0f)
	{
		num.assign(3, 0);
		PRECISION = 1000.0f;
	}
	else if(worldscale < 5.0f)
	{
		num.assign(2, 0);
		PRECISION = 100.f;
	}
	
	bool neg = false;
	if(number < 0)
	{
		number *= -1;
		neg = true;
	}
	int tmp1 = (int) number;
	if(number != (int) number)
	{
		int tmp2 = ((int) (number * PRECISION)) - (((int) number) * PRECISION), cnt = 0;
		while(tmp2 > 9)
		{
			num[cnt++] = tmp2 % 10;
			tmp2 /= 10;
		}
		num[cnt++] = tmp2;
	}
	num.push_back(-1);
	while(tmp1 > 9)
	{
		num.push_back(tmp1 % 10);
		tmp1 /= 10;
	}
	num.push_back(tmp1);
	bool coda = true, decimal = true;
	std::vector<double> v;
	for(int i = 0; i < num.size(); ++i)
	{
		if(num[i] == -1)
		{
			decimal = false;
			if(v.size() == 0) continue;
			for(int j = 0; j < nums[10].size(); ++j)
			{
				int el = nums[10][j] * scale;
				if(j % 2 == 0) el += scale * (num.size() - i - 1) + scale * 0.5;
				v.push_back(el);
			}
			continue;
		}
		if(num[i] == 0 && coda && decimal) continue;
		for(int j = 0; j < nums[num[i]].size(); ++j)
		{
			int el = nums[num[i]][j] * scale;
			if(j % 2 == 0) el += scale * (num.size() - i - 1) + ((decimal) ? 0 : scale * 0.5);
			v.push_back(el);
		}
		coda = false;
	}
	if(neg)
	{
		for(int j = 0; j < nums[11].size(); ++j)
		{
			int el = nums[11][j] * scale;
			if(j % 2 == 0) el -= scale / 2;
			v.push_back(el);
		}
	}
	if(option == ALIGN_RIGHT)
	{
		double M = 0;
		for(int i = 0; i < v.size(); i += 2) M = std::max(M, v[i]);
		for(int i = 0; i < v.size(); i += 2) v[i] -= M * 1.5f;
		for(int i = 1; i < v.size(); i += 2) v[i] -= scale / 2;
	}
	else if(option == CENTERED)
	{
		double M = 0;
		for(int i = 0; i < v.size(); i += 2) M = std::max(M, v[i]);
		for(int i = 0; i < v.size(); i += 2) v[i] -= M / 2;
		for(int i = 1; i < v.size(); i += 2) v[i] -= scale * 1.5f;
	}
	return v;
}

void Init(int width, int height, const std::string& title)
{
	glfwInit();
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	window = glfwCreateWindow(width, height, title.data(), NULL, NULL);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);
	glewInit();
	ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    
	glfwSetCursorPosCallback(window, MousePosProc);
	glfwSetScrollCallback(window, MouseScrollProc);
	glfwSetMouseButtonCallback(window, MouseButtonProc);
//	glfwSetKeyCallback(window, KeyboardProc);
	glDebugMessageCallback(ErrorProc, 0);
	
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT,  GL_NICEST);
	
    glGenBuffers(1, &vertexbuffer);
    glGenBuffers(1, &axisbuffer);
    glGenBuffers(1, &numbersbuffer);
    glGenBuffers(1, &deepgridbuffer);
    glGenBuffers(1, &gridbuffer);
    glGenBuffers(1, &integratebuffer);
    glGenBuffers(1, &guidelinesbuffer);
    glGenBuffers(1, &parabolasbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void Destroy()
{
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteBuffers(1, &axisbuffer);
	glDeleteBuffers(1, &numbersbuffer);
	glDeleteBuffers(1, &gridbuffer);
	glDeleteBuffers(1, &deepgridbuffer);
	glDeleteBuffers(1, &integratebuffer);
	glDeleteBuffers(1, &guidelinesbuffer);
	glDeleteBuffers(1, &parabolasbuffer);
	glfwTerminate();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void APIENTRY ErrorProc(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam)
{
	std::cout << "Error : " << message << "\n";
	if(severity == GL_DEBUG_SEVERITY_HIGH)
	{
		std::cout << "ABORTING ...";
		std::abort();
	}
}

void KeyboardProc(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	
}

void MousePosProc(GLFWwindow* window, double x, double y)
{
	currMousePos = {x, y};
	if(drag)
	{
		xoffset += (lastMousePos.first - currMousePos.first) * scale;
		yoffset += (currMousePos.second - lastMousePos.second) * scale;
		lastMousePos = currMousePos;
	}
}

void MouseScrollProc(GLFWwindow* window, double xoff, double yoff)
{
	double new_scale = scale * (1.0 + yoff / 10.0f);
	if(new_scale < 0.001f)
		new_scale = 0.001f;
	
	if(precision_scale && new_scale <= 1.0f)
		scale = new_scale;
	else if(!precision_scale && new_scale >= 1.0f)
		scale = new_scale;
	else
		scale = 1.0f;
}

void MouseButtonProc(GLFWwindow* window, int button, int action, int mods)
{
	if(button == GLFW_MOUSE_BUTTON_MIDDLE)
	{
		if(action == GLFW_PRESS)
		{
			drag = true;
			lastMousePos = currMousePos;
		}
		else if(action == GLFW_RELEASE)
		{
			drag = false;
		}
	}
}

