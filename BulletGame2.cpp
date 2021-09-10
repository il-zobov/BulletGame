#include <vector>
#include <list> 
#include <algorithm>
#include <cmath>
#include <limits>

#include <chrono>
#include <thread>
#include <mutex>

using namespace std;

struct float2
{
	float x;
	float y;
};

struct LineSegment
{
	float2 first;
	float2 second;
};

struct Box
{
	float2 first;
	float2 second;
};

/* smallest such that 1.0+FLT_EPSILON != 1.0 */
bool is_equal( float x, float y)
{
	return fabs(x - y) < std::numeric_limits<float>::epsilon(); 
}

static double crossProduct(const float2 &a,const float2 &b) 
{
	return a.x * b.y - b.x * a.y;
}

static bool doBoundingBoxesIntersect(const Box &a,const Box &b) 
{
	// may be loose of data when == happen
	return a.first.x <= b.second.x && a.second.x >= b.first.x && a.first.y <= b.second.y && a.second.y >= b.first.y;
}

static bool isPointOnLine(const LineSegment &a,const float2 &b)
{
	float2 zeroPoint = { 0, 0 };
	LineSegment aTmp = { zeroPoint, { a.second.x - a.first.x, a.second.y - a.first.y } };
	float2 bTmp = { b.x - a.first.x, b.y - a.first.y };
	double r = crossProduct(aTmp.second, bTmp);
	return abs(r) <  std::numeric_limits<float>::epsilon();
}

static bool isPointRightOfLine(const LineSegment &a,const float2 &b) 
{
	float2 zeroPoint = { 0, 0 };
	LineSegment aTmp = { zeroPoint, { a.second.x - a.first.x, a.second.y - a.first.y } };
	float2 bTmp = { b.x - a.first.x, b.y - a.first.y };
	return crossProduct(aTmp.second, bTmp) < 0;
}

static bool lineSegmentTouchesOrCrossesLine (const LineSegment &a,const LineSegment &b)
{
	return isPointOnLine(a, b.first)
		|| isPointOnLine(a, b.second)
		|| (isPointRightOfLine(a, b.first) ^ isPointRightOfLine(a,
		b.second));
}

Box getBoundingBox(const LineSegment &ls)
{
	Box result;
	result.first = { std::min(ls.first.x, ls.second.x), std::min(ls.first.y, ls.second.y) };
	result.second = { std::max(ls.first.x, ls.second.x), std::max(ls.first.y, ls.second.y) };
	return result;
}

// Check if to line segments intersect 
static bool doLinesIntersect(const LineSegment &a,const LineSegment &b) {
	Box box1 = getBoundingBox(a);
	Box box2 = getBoundingBox(b);
	return doBoundingBoxesIntersect(box1, box2)
		&& lineSegmentTouchesOrCrossesLine(a, b)
		&& lineSegmentTouchesOrCrossesLine(b, a);
}

/** You know that lines a and b have an intersection and now you     want to get it! */
// if the intersection is a line then return is line segment 
// otherwise x1=x2,y1=y2, so it is a point
LineSegment getIntersection(LineSegment a,LineSegment b) {
	/* the intersection [(x1,y1), (x2, y2)]
	it might be a line or a single point. If it is a point,        then x1 = x2 and y1 = y2.  */
	float x1, y1, x2, y2;
	if (is_equal(a.first.x, a.second.x)) {
		x1 = a.first.x;
		x2 = x1;
		if (is_equal(b.first.x, b.second.x))
		{
			if (a.first.y > a.second.y)
			{
				a = { a.second, a.first };
			}
			if (b.first.y > b.second.y)
			{
				b = { b.second, b.first };
			}
			if (a.first.y > b.first.y)
			{
				LineSegment temp = a;
				a = b;
				b = temp;
			}

			y1 = b.first.y;
			y2 = std::min(a.second.y, b.second.y);
		}
		else
		{
			float m, t;
			m = (b.first.y - b.second.y) / (b.first.x - b.second.x);
			t = b.first.y - m*b.first.x;
			y1 = m*x1 + t;
			y2 = y1;

		}
	}
	else if (is_equal(b.first.x, b.second.x))
	{
		x1 = b.first.x;
		x2 = x1;

		LineSegment temp = a;
		a = b;
		b = temp;
		float m, t;
		m = (b.first.y - b.second.y) / (b.first.x - b.second.x);
		t = b.first.y - m*b.first.x;
		y1 = m*x1 + t;
		y2 = y1;
	}
	else
	{
		float ma, mb, ta, tb;
		ma = (a.first.y - a.second.y) / (a.first.x - a.second.x);
		mb = (b.first.y - a.second.y) / (b.first.x - b.second.x);
		ta = a.first.y - ma*a.first.x;
		tb = b.first.y - mb*b.first.x;
		if (is_equal(ma, mb))
		{
			if (a.first.x > a.second.x)
			{
				a = { a.second, a.first };
			}
			if (b.first.x > b.second.x)
			{
				b = { b.second, b.first };
			}
			if (a.first.x > b.first.x)
			{
				LineSegment temp = a;
				a = b;
				b = temp;
			}
			x1 = b.first.x;
			x2 = min(a.second.x, b.second.x);
			y1 = ma*x1 + ta;
			y2 = ma*x2 + ta;
		}
		else
		{
			x1 = (tb - ta) / (ma - mb);
			y1 = ma*x1 + ta;
			x2 = x1;
			y2 = y1;
		}
	}
	return{ { x1, y1 }, { x2, y2 } };
}

float2 reflectVector(const float2 *pa, const float2  *pb, const float2 *ps)
{
	float2 n;
	float2 result;
	n.x = (pa->y - pb->y);
	n.y = (pb->x - pa->x);
	float len = sqrt(n.x*n.x + n.y*n.y);
	n.x /= len;
	n.y /= len;
	float dot2 = 2 * (n.x*ps->x + n.y*ps->y);
	result.x = ps->x - dot2*n.x;
	result.y = ps->y - dot2*n.y;
	return result;
}

class Bullet
{
public:
	// time of bulet start flight
	float fireTime;
	float lifeTime;
	// time of larst bullet position updated 
	float updateTime;
	// distance betwen the bullet and the wall thet the bullet hits
	float speed;
	float2 poz;
	float2 dir;

	Bullet(float startTime, float life, float bulSpeed, float2 pozition, float2 direction);
	~Bullet();
	LineSegment getCurrenBulletFlyLine(const float &currTime);

private:

	float2 normalizeVector( const float2 &vector);
};

Bullet::Bullet(float startTime, float life, float bulSpeed, float2 pozition, float2 direction)
{	
	fireTime = startTime;
	lifeTime = life;
	speed = bulSpeed;
	poz = pozition;
	dir = direction;
	updateTime = startTime;
	
	printf("Bullet Started fligtNomber %.f\n", fireTime);
}

Bullet::~Bullet()
{
	 printf("Bullet dead fligtNomber %.f\n", fireTime);
}

LineSegment Bullet::getCurrenBulletFlyLine(const float &deltaTime)
{
	float flighTime = deltaTime;
	float2 normDir = normalizeVector(dir);
	float xcoord = poz.x + speed*normDir.x;
	float  ycoord = poz.y + speed*normDir.y;
	return{ poz, { xcoord, ycoord } };
}

float2 Bullet::normalizeVector(const float2 &vector)
{
	float2 result = { 0, 0 };
	float modul = sqrt(pow(vector.x, 2) + pow(vector.y, 2));
	if (!is_equal(modul, 0))
		result = { vector.x / modul, vector.y / modul };
	return result;
}

bool sortByTime(const Bullet* lhs, const Bullet* rhs) {
	return (lhs->fireTime) < (rhs->fireTime);
}

// each cross is a possible intersection of wall and bullet
struct Cross
{
	Bullet* bul;
	LineSegment* wall;
	float timeTillIntersect;
	float2 intersectPoint;
};

bool sortByPointer(const Cross &lhs, const Cross &rhs) {
	return  (lhs.wall) < (rhs.wall);
}

bool uniq(const Cross &lhs, const Cross &rhs) {
	return lhs.bul ==rhs.bul;
// could be a missing of points here
}

bool linesEquel( LineSegment a, LineSegment b)
{
	return (is_equal(a.first.x, b.first.x) && is_equal(a.first.y, b.first.y) && is_equal(a.second.x, b.second.x) && is_equal(a.second.y, b.second.y));
}

class BulletManager
{
public:
	vector<Bullet*> incoming;
	vector<Bullet*> inFlight;
	vector<Cross> intersections;
	list<LineSegment> walls;

	BulletManager(list<LineSegment> wall);

	// puting the bullets to incoming collection
	void fire(float2 pos, float2 dir, float speed, float time, float life_time, mutex &mtx);
	void BulletManager::update(float time, mutex &mtx);

private:
	void printBulVect(const vector<Bullet*> &vector)
	{
		for (auto it = vector.begin(); it != vector.end(); ++it)
		{
			printf("Vector  %.f\n", (**it).fireTime);
		}
	}

	float distenceBetwPoints(float2 &a, float2 &b)
	{
		float ac = b.x - a.x;
		float bc = b.y - a.y;
		return sqrt((ac*ac) + (bc*bc));
	}
};

BulletManager::BulletManager(list<LineSegment> wall)
{
	walls = wall;
}

// puting the bullets to incoming collection
void BulletManager::fire(float2 pos, float2 dir, float speed, float time, float life_time, mutex &mtx)
{
	mtx.lock();
		incoming.push_back( new Bullet(time, life_time, speed, pos, dir));
	mtx.unlock();
	this_thread::sleep_for(std::chrono::milliseconds(50));
}

//  put bullet to flight
// who crossed the wall
// who is dead 
void BulletManager::update(float time, mutex &mtx)
{
	mtx.lock();
		// sorting vector of incoming bullets
		sort(incoming.begin(), incoming.end(), sortByTime);
		printBulVect(incoming);
		for (auto it = incoming.begin(); it != incoming.end(); it++)
		{
			if ((**it).fireTime > time)
			{
				inFlight.insert(inFlight.end(), incoming.begin(), it);
				incoming.erase(incoming.begin(), it);
				break;
			}
		}
	mtx.unlock();

	for (auto wallIt = walls.begin(); wallIt != walls.end(); ++wallIt){
		for (auto bullIt = inFlight.begin(); bullIt != inFlight.end(); ++bullIt)
		{
			LineSegment flyline = (*bullIt)->getCurrenBulletFlyLine(time - (*bullIt)->updateTime);
			// making a vector of possible intersections
			if (doLinesIntersect(flyline,*wallIt))
			{
				float2 inersectionPoint = getIntersection(flyline, *wallIt).first;
				float intersectionTime = (distenceBetwPoints(inersectionPoint, (*bullIt)->poz)) / (*bullIt)->speed;
				// check the bullet is alive when it hits the wall

				// wall is already in a vector
				bool wallPresent = false;
				for (auto it = intersections.begin(); it != intersections.end(); ++it)
				{
					// if wall present then check who is the firs bullet to hit it
					if (it->wall == &(*wallIt))
					{
						wallPresent = true;
						if (it->timeTillIntersect > intersectionTime)
						{
							it->bul = *bullIt;
							it->timeTillIntersect = intersectionTime; 
							it->intersectPoint = inersectionPoint;
						}
					}
				}
				// if no wall present in vector then, thet is the first bullet in it
				if (!wallPresent)
				{
					intersections.push_back({ *bullIt, &(*wallIt), intersectionTime, inersectionPoint });
				}
			}
			//else..
		}
	}
	// sort vector for excludin similar instances
	sort(intersections.begin(), intersections.end(), sortByPointer);
	// exclude same bullets
	intersections.resize(unique(intersections.begin(), intersections.end(), uniq) - intersections.begin());

	for (auto it = intersections.begin(); it != intersections.end(); ++it)
	{
		// chek if the bulet time till intersect is ok
		// reflect vector
		// find coordinat 
		// earse wall 
		// check if bullet still exists at the end

		float timeOfBulletDeath= (*it->bul).fireTime + (*it->bul).lifeTime;
		(*it->bul).dir = reflectVector(&(*it->wall).first, &(*it->wall).second, &(*it->bul).dir);
		if (timeOfBulletDeath > time)
		{
			// cross wall and bullet
			(*it->bul).poz = it->intersectPoint;
			// now the bullet reflects and fly to until update time ends, so there is new pozition
			(*it->bul).poz = (*it->bul).getCurrenBulletFlyLine((time - ((*it->bul).updateTime) + it->timeTillIntersect)).second;
			(*it->bul).updateTime = time;
			
			// remove this wall from list
			walls.remove_if(bind2nd(ptr_fun(linesEquel), (*it->wall)));
		}
		//else...
	}
	// clear the vector of intersections
	intersections.clear();

	for (auto bullIt = inFlight.begin(); bullIt != inFlight.end(); ++bullIt)
	{
		if (((*bullIt)->lifeTime + (*bullIt)->fireTime)< time)
		{
			delete	*bullIt;
			*bullIt = NULL;
		}
	}
	inFlight.erase(remove(inFlight.begin(), inFlight.end(), (Bullet*)NULL), inFlight.end());
}

int main()
{
	mutex mtx;

	list<LineSegment> wall = {
		{ { 5, 0 }, { 5, 10 } },
		{ { 1, 0 }, { 1, 10 } },
		{ { 15, 0 }, { 15, 10 } },
		{ { 15, 0 }, { 15, 10 } }
	};



	BulletManager manager = BulletManager(wall);
	manager.fire({ 0, 0 }, { 1, 1 }, 60, 1, 30,mtx);
	manager.fire({ 0, 0 }, { 1, 1 }, 60, 6, 40, mtx);
	manager.fire({ 0, 0 }, { 1, 1 }, 50, 5, 50, mtx);
	manager.fire({ 0, 0 }, { 1, 1 }, 120, 2, 60, mtx);

	
	manager.update(5, mtx);



	getchar();
	return 0;
}

