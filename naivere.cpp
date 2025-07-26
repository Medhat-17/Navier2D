#include <vector>
#include <cmath>
#include <unordered_map>
#include <iostream>

struct Vector2
{
    float x, y;
    Vector2(float x = 0.0f, float y = 0.0f) : x(x), y(y) {}
    Vector2 operator+(const Vector2 &other) const { return Vector2(x + other.x, y + other.y); }
    Vector2 operator-(const Vector2 &other) const { return Vector2(x - other.x, y - other.y); }
    Vector2 operator*(float scalar) const { return Vector2(x * scalar, y * scalar); }
    Vector2 operator/(float scalar) const { return Vector2(x / scalar, y / scalar); }
    Vector2 &operator+=(const Vector2 &other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }
    float dot(const Vector2 &other) const { return x * other.x + y * other.y; }
    float magnitude() const { return std::sqrt(x * x + y * y); }
    float magnitudeSquared() const { return x * x + y * y; }
};

// Particle class for SPH
struct Particle
{
    Vector2 position;
    Vector2 velocity;
    Vector2 force;
    float density;
    float pressure;
    Particle(Vector2 pos) : position(pos), velocity(0, 0), force(0, 0), density(0), pressure(0) {}
};

// SpatialHash for efficient neighbor search
class SpatialHash
{
private:
    std::unordered_map<int, std::vector<int>> grid;
    float cellSize;

    int getCellKey(float x, float y) const
    {
        int i = static_cast<int>(x / cellSize);
        int j = static_cast<int>(y / cellSize);
        return i + j * 10000; // Simple hash for 2D grid
    }

public:
    SpatialHash(float size) : cellSize(size) {}

    void clear() { grid.clear(); }

    void insert(int particleIdx, const Vector2 &pos)
    {
        int key = getCellKey(pos.x, pos.y);
        grid[key].push_back(particleIdx);
    }

    std::vector<int> queryNeighbors(const Vector2 &pos, float radius) const
    {
        std::vector<int> neighbors;
        int i = static_cast<int>(pos.x / cellSize);
        int j = static_cast<int>(pos.y / cellSize);
        for (int di = -1; di <= 1; ++di)
        {
            for (int dj = -1; dj <= 1; ++dj)
            {
                int key = (i + di) + (j + dj) * 10000;
                auto it = grid.find(key);
                if (it != grid.end())
                {
                    neighbors.insert(neighbors.end(), it->second.begin(), it->second.end());
                }
            }
        }
        return neighbors;
    }
};

// SPHSolver class for fluid simulation
class SPHSolver
{
private:
    std::vector<Particle> particles;
    SpatialHash spatialHash;
    float h;                      // Smoothing radius
    float h2;                     // Smoothing radius squared
    float mass;                   // Particle mass
    float restDensity;            // Reference density
    float viscosity;              // Viscosity coefficient
    float k;                      // Pressure constant (stiffness)
    float dt;                     // Time step
    Vector2 boundsMin, boundsMax; // Simulation boundaries

    // SPH kernels
    float poly6Kernel(float r2) const
    {
        if (r2 > h2)
            return 0.0f;
        float term = h2 - r2;
        return (315.0f / (64.0f * M_PI * std::pow(h, 9))) * term * term * term;
    }

    Vector2 spikyKernelGradient(Vector2 r, float rLen) const
    {
        if (rLen > h || rLen < 0.0001f)
            return Vector2(0, 0);
        float term = h - rLen;
        float coeff = -45.0f / (M_PI * std::pow(h, 6)) * term * term;
        return r * (coeff / rLen);
    }

    float viscosityKernelLaplacian(float rLen) const
    {
        if (rLen > h)
            return 0.0f;
        return (45.0f / (M_PI * std::pow(h, 6))) * (h - rLen);
    }

    // Compute density and pressure
    void computeDensityPressure()
    {
        for (size_t i = 0; i < particles.size(); ++i)
        {
            particles[i].density = 0.0f;
            auto neighbors = spatialHash.queryNeighbors(particles[i].position, h);
            for (int j : neighbors)
            {
                Vector2 r = particles[i].position - particles[j].position;
                float r2 = r.magnitudeSquared();
                particles[i].density += mass * poly6Kernel(r2);
            }
            // Tait equation for pressure
            particles[i].pressure = k * (std::pow(particles[i].density / restDensity, 7) - 1.0f);
        }
    }

    // Compute forces (pressure, viscosity, external)
    void computeForces()
    {
        for (size_t i = 0; i < particles.size(); ++i)
        {
            particles[i].force = Vector2(0, 0);
            auto neighbors = spatialHash.queryNeighbors(particles[i].position, h);
            Vector2 pressureForce(0, 0);
            Vector2 viscosityForce(0, 0);

            for (int j : neighbors)
            {
                if (i == j)
                    continue;
                Vector2 r = particles[i].position - particles[j].position;
                float rLen = r.magnitude();
                if (rLen < h && rLen > 0.0001f)
                {
                    // Pressure force
                    float pTerm = -mass * (particles[i].pressure + particles[j].pressure) / (2.0f * particles[j].density);
                    pressureForce += spikyKernelGradient(r, rLen) * pTerm;
                    // Viscosity force
                    Vector2 vDiff = particles[j].velocity - particles[i].velocity;
                    viscosityForce += vDiff * (mass * viscosity * viscosityKernelLaplacian(rLen) / particles[j].density);
                }
            }
            particles[i].force += pressureForce + viscosityForce;
            // External force (gravity)
            particles[i].force += Vector2(0, -9.81f * mass);
        }
    }

    // Handle boundary collisions
    void handleBoundaries()
    {
        for (auto &p : particles)
        {
            float restitution = 0.9f;
            if (p.position.x < boundsMin.x)
            {
                p.position.x = boundsMin.x;
                p.velocity.x = -p.velocity.x * restitution;
            }
            if (p.position.x > boundsMax.x)
            {
                p.position.x = boundsMax.x;
                p.velocity.x = -p.velocity.x * restitution;
            }
            if (p.position.y < boundsMin.y)
            {
                p.position.y = boundsMin.y;
                p.velocity.y = -p.velocity.y * restitution;
            }
            if (p.position.y > boundsMax.y)
            {
                p.position.y = boundsMax.y;
                p.velocity.y = -p.velocity.y * restitution;
            }
        }
    }

public:
    SPHSolver(float smoothingRadius, float particleMass, float restDens, float visc, float pressureK, float timeStep)
        : spatialHash(smoothingRadius), h(smoothingRadius), h2(h * h), mass(particleMass),
          restDensity(restDens), viscosity(visc), k(pressureK), dt(timeStep),
          boundsMin(0, 0), boundsMax(10, 10) {}

    void addParticle(Vector2 pos)
    {
        particles.emplace_back(pos);
        spatialHash.insert(particles.size() - 1, pos);
    }

    void setBounds(Vector2 min, Vector2 max)
    {
        boundsMin = min;
        boundsMax = max;
    }

    void update()
    {
        // Update spatial hash
        spatialHash.clear();
        for (size_t i = 0; i < particles.size(); ++i)
        {
            spatialHash.insert(i, particles[i].position);
        }
        // Compute densities and pressure
        computeDensityPressure();
        // Compute forces
        computeForces();
        // Update velocities and positions (Leapfrog integration)
        for (auto &p : particles)
        {
            p.velocity += p.force * (dt / mass);
            p.position += p.velocity * dt;
        }
        // Handle boundaries
        handleBoundaries();
    }

    // For debugging/output
    void printState() const
    {
        for (size_t i = 0; i < particles.size(); ++i)
        {
            std::cout << "Particle " << i << ": Pos(" << particles[i].position.x << ", "
                      << particles[i].position.y << "), Vel(" << particles[i].velocity.x << ", "
                      << particles[i].velocity.y << "), Density(" << particles[i].density << ")\n";
        }
    }
};

// Main function to test the simulation
int main()
{
    // Initialize SPH solver: h=0.5, mass=0.02, restDensity=1000, viscosity=0.1, k=1000, dt=0.01
    SPHSolver solver(0.5f, 0.02f, 1000.0f, 0.1f, 1000.0f, 0.01f);

    // Add particles in a 5x5 grid
    for (float x = 2.0f; x <= 4.0f; x += 0.2f)
    {
        for (float y = 2.0f; y <= 4.0f; y += 0.2f)
        {
            solver.addParticle(Vector2(x, y));
        }
    }

    // Set simulation bounds
    solver.setBounds(Vector2(0, 0), Vector2(10, 10));

    // Simulate for 100 steps
    for (int i = 0; i < 100; ++i)
    {
        solver.update();
        solver.printState();
    }

    return 0;
}