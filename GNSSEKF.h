#pragma once
#define DeltaT 1.0

struct GNSSEKF
{
	Matrix * X;//X Vx Ax Y Vy Ay Z Vz Az T dT
	Matrix * Xp;
	Matrix * Dx;
	Matrix * Dp;
	Matrix * F;
	Matrix * T;

	Matrix * Ft;
	Matrix * Tt;
	Matrix * De;

	Matrix * Z;
	Matrix * H;
	Matrix * Dz;
	Matrix * Zp;

	Matrix * K;
	Matrix * V;
};
GNSSEKF EKFCreate();
void EKFProcess(
	GNSSEKF * ekf,
	double * fake_distance,
	double * fake_distance_err,
	SpaceLocation * satellites_position,
	int sat_num
);

