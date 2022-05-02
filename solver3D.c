#define IX(i,j,k) ((i)+(N+2)*(j)+(N+2)*(N+2)*(k))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) { for ( k=1 ; k<=N ; k++ ){
#define END_FOR }}}

void add_source(int N, float* x, float* s, float dt)
{
	int i, size = (N + 2) * (N + 2)*(N+2);
	for (i = 0; i < size; i++) x[i] += dt * s[i];
}

void set_bnd(int N, int b, float* x)
{
	int i,z;
	for (z = 1; z <= N; z++) {
		for (i = 1; i <= N; i++) {
			x[IX(0, i, z)] = b == 1 ? -x[IX(1, i, z)] : x[IX(1, i, z)];
			x[IX(N + 1, i, z)] = b == 1 ? -x[IX(N, i, z)] : x[IX(N, i, z)];
			x[IX(i, 0, z)] = b == 2 ? -x[IX(i, 1, z)] : x[IX(i, 1, z)];
			x[IX(i, N + 1, z)] = b == 2 ? -x[IX(i, N, z)] : x[IX(i, N, z)];
			x[IX(i, z, 0)] = b == 3 ? -x[IX(i, z, 1)] : x[IX(i, z, 1)];
			x[IX(i, z, N + 1)] = b == 3 ? -x[IX(i, z, N)] : x[IX(i, z, N)];
		}
	}

	x[IX(0, 0, 0)] = 1.0f/3 * (x[IX(1, 0,0)] + x[IX(0, 1,0)] + x[IX(0, 0, 1)]);
	x[IX(N+1, 0, 0)] = 1.0f / 3 * (x[IX(N, 0,0)] + x[IX(N+1,1,0)] + x[IX(N + 1, 0, 1)]);
	x[IX(0, N + 1,0)] = 1.0f / 3 * (x[IX(0, N, 0)] + x[IX(1, N+1, 0)] + x[IX(0, N + 1, 1)]);
	x[IX(0, 0, N + 1)] = 1.0f / 3 * (x[IX(0, 0, N)] + x[IX(1, 0, N+1)] + x[IX(0, 1, N + 1)]);
	x[IX(N+1, N+1, 0)] = 1.0f / 3 * (x[IX(N, N+1, 0)] + x[IX(N+1, N, 0)] + x[IX( N + 1, N+1, 1)]);
	x[IX(N+1, 0, N + 1)] = 1.0f / 3 * (x[IX(N, 0, N+1)] + x[IX(N+1, 1, N + 1)] + x[IX(N+1, 0, N)]);
	x[IX(0, N+1, N + 1)] = 1.0f / 3 * (x[IX(1, N+1, N+1)] + x[IX(0, N, N + 1)] + x[IX(0, N+1, N)]);
	x[IX(N + 1, N + 1, N+1)] = 1.0f / 3 * (x[IX(N, N + 1, N+1)] + x[IX(N + 1, N, N+1)]+ x[IX(N + 1, N+1, N)]);
}

void lin_solve(int N, int b, float* x, float* x0, float a, float c)
{
	int i, j, m,k;

	for (m = 0; m < 20; m++) {
		FOR_EACH_CELL
			x[IX(i, j,k)] = (x0[IX(i, j,k)] + a * (x[IX(i - 1, j,k)] + x[IX(i + 1, j,k)] + x[IX(i, j - 1,k)] + x[IX(i, j + 1,k)+ x[IX(i, j, k+1)] + x[IX(i, j, k-1)])) / c;
		END_FOR
			set_bnd(N, b, x);
	}
}

void diffuse(int N, int b, float* x, float* x0, float diff, float dt)
{
	float a = dt * diff * N * N*N;
	lin_solve(N, b, x, x0, a, 1 + 6 * a);
}

void advect(int N, int b, float* d, float* d0, float* u, float* v, float* w, float dt)
{
	int i, j, k,  i1, j1, k1, i2,j2,k2;
	float x, y, z,dt0;

	dt0 = dt * N;
	FOR_EACH_CELL
		x = i - dt0 * u[IX(i, j,k)];
		y = j - dt0 * v[IX(i, j,k)];
		z = k - dt0 * w[IX(i, j, k)];

		if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f;
		i1 = (int)x; i2 = i1 + 1;
		if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; 
		j1 = (int)y; j2 = j1 + 1;
		if (z < 0.5f) z = 0.5f; if (z > N + 0.5f) z = N + 0.5f;
		k1 = (int)y; k2 = k1 + 1; 

		d[IX(i, j, k)] = d0[IX(i1, j1, k1)] * (i2 - x) * (j2 - y) * (k2 - z) + d0[IX(i2, j1, k1)] * (x - i1) * (j2 - y) * (k2 - z)
			+ d0[IX(i2, j1, k2)] * (x - i1) * (j2 - y) * (z - k1) + d0[IX(i2, j2, k1)] * (x - i1) * (y - j1) * (k2 - z)
			+ d0[IX(i2, j2, k2)] * (x - i1) * (y - j1) * (z - k1);
	END_FOR
	set_bnd(N, b, d);
}

void project(int N, float* u, float* v, float* w, float* p, float* div)
{
	int i, j,k;

	FOR_EACH_CELL
		div[IX(i, j,k)] = -0.5f * (u[IX(i + 1, j,k)] - u[IX(i - 1, j,k)] + v[IX(i, j + 1,k)] - v[IX(i, j - 1,k)]
			+ w[IX(i, j, k+1)] - w[IX(i, j, k+1)]) / N;
	p[IX(i, j,k)] = 0;
	END_FOR
	set_bnd(N, 0, div); set_bnd(N, 0, p);

	lin_solve(N, 0, p, div, 1, 6);

	FOR_EACH_CELL
		u[IX(i, j,k)] -= 0.5f * N * (p[IX(i + 1, j,k)] - p[IX(i - 1, j,k)]);
		v[IX(i, j,k)] -= 0.5f * N * (p[IX(i, j + 1,k)] - p[IX(i, j - 1,k)]);
		w[IX(i, j, k)] -= 0.5f * N * (p[IX(i, j, k+1)] - p[IX(i, j, k-1)]);
	END_FOR
		set_bnd(N, 1, u); set_bnd(N, 2, v);
}

void dens_step(int N, float* x, float* x0, float* u, float* v, float* w,float diff, float dt)
{
	add_source(N, x, x0, dt);
	SWAP(x0, x); diffuse(N, 0, x, x0, diff, dt);
	SWAP(x0, x); advect(N, 0, x, x0, u, v, w, dt);
}

void vel_step(int N, float* u, float* v, float* w, float* u0, float* v0, float* w0, float visc, float dt)
{
	add_source(N, u, u0, dt); add_source(N, v, v0, dt); add_source(N, w, w0, dt);
	SWAP(u0, u); diffuse(N, 1, u, u0, visc, dt);
	SWAP(v0, v); diffuse(N, 2, v, v0, visc, dt);
	SWAP(w0, w); diffuse(N, 2, w, w0, visc, dt);
	project(N, u, v, w, u0, v0);
	SWAP(u0, u); SWAP(v0, v); SWAP(w0, w);
	advect(N, 1, u, u0, u0, v0, w0, dt); advect(N, 2, v, v0, u0, v0, w0,dt); advect(N, 2, w, v0, u0, v0, w0, dt);
	project(N, u, v,w, u0, v0);
}

