using System;
using System.Linq;
using UnityEngine;

namespace LBMNamespace
{
    public class LBM
    {
        // Private fixed parameters
        private const double p0 = 0;
        private const double nu_LB = 0.01;
        private const double cs2 = 1.0 / 3.0;
        private const double tau_LB = (nu_LB + 0.5) / cs2;

        private int nx, ny, currentIteration, numberOfIterations;
        private double uMax, dx, dt;
        public double[] density;
        public double[,] velocity;
        private int[] WallMap;
        private double[,] df1, df2;
        private double omega;
        private double[] weights = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
        private double[,] velocitySet = {
            { 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0 },
            { 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0 }
        };

        private const int flowID = 0;
        private const int wallID = 1;
        private const int inflowID = 2;
        private const int outflowID = 3;

        // Define Cx, Cy, and R for the MapInit method
        private const double Cx = 1.0;
        private const double Cy = 1.0;
        private const double R = 0.5;

        public LBM(double uMax, double nu_phys, double Lx, double Ly, double T, double dx)
        {
            Debug.Log("LBM constructor");
            this.uMax = uMax;
            this.dx = dx;

            nx = (int)Math.Ceiling(Lx / dx);
            ny = (int)Math.Ceiling(Ly / dx);

            dt = nu_LB / nu_phys * dx * dx;
            currentIteration = 0;
            numberOfIterations = (int)Math.Ceiling(T / dt);

            density = new double[nx * ny];
            velocity = new double[nx * ny, 2];
            WallMap = new int[nx * ny];

            df1 = new double[nx * ny, 9];
            df2 = new double[nx * ny, 9];

            omega = cs2 / (nu_LB + 0.5);
        }

        private double[] Equilibrium(double i_rho, double i_u, double i_v)
        {
            Debug.Log("Equilibrium");
            double[] eq = new double[9];
            for (int idx = 0; idx < 9; idx++)
            {
                eq[idx] = i_rho * weights[idx] * (1.0 + 3.0 * (velocitySet[0, idx] * i_u + velocitySet[1, idx] * i_v)
                    + 4.5 * Math.Pow(velocitySet[0, idx] * i_u + velocitySet[1, idx] * i_v, 2)
                    - 1.5 * (i_u * i_u + i_v * i_v));
            }
            return eq;
        }

        private (double, double[]) ComputeMacro(double[] i_df)
        {
            Debug.Log("ComputeMacro");
            double rho_loc = i_df.Sum();
            double vx_loc = i_df.Select((df, idx) => df * velocitySet[0, idx]).Sum() / rho_loc;
            double vy_loc = i_df.Select((df, idx) => df * velocitySet[1, idx]).Sum() / rho_loc;
            return (rho_loc, new double[] { vx_loc, vy_loc });
        }

        private (double[], double, double[]) InflowBoundary(int i_x, int i_y, double[] i_df)
        {
            Debug.Log("InflowBoundary");
            var (xf, yf) = PositionCellCenter(i_x, i_y);
            double[] loc_df_plus = Streaming(i_df, i_x + 1, i_y);
            var (rho_plus, _) = ComputeMacro(loc_df_plus);
            double no1Ocs2 = 3.0;
            double InputFactor = no1Ocs2 * Phys2LBMVelocity(Phys2LBMVelocity(nu_LB * 8.0 * uMax / Math.Pow(ny * dx, 2)));
            var (phys_vx, phys_vy) = InflowVelocity(yf);
            double loc_vx = Phys2LBMVelocity(phys_vx);
            double loc_vy = Phys2LBMVelocity(phys_vy);
            double loc_rho = rho_plus - InputFactor;
            return (Equilibrium(loc_rho, loc_vx, loc_vy), loc_rho, new double[] { loc_vx, loc_vy });
        }

        private double[] Bounceback(double[] i_df)
        {
            Debug.Log("Bounceback");
            return new double[] { i_df[0], i_df[3], i_df[4], i_df[1], i_df[2], i_df[7], i_df[8], i_df[5], i_df[6] };
        }

        private double[] OutflowCondition(int i_x, int i_y, double[] i_df)
        {
            Debug.Log("OutflowCondition");
            return i_df;
        }

        private double[] Streaming(double[] i_df, int i_x, int i_y)
        {
            Debug.Log("Streaming");
            double[] out_df = new double[9];
            out_df[0] = i_df[D1Idx(i_x, i_y)];
            int xm = Math.Max(0, i_x - 1);
            int ym = Math.Max(0, i_y - 1);
            int xp = Math.Min(nx - 1, i_x + 1);
            int yp = Math.Min(ny - 1, i_y + 1);
            out_df[1] = i_df[D1Idx(xm, i_y)];
            out_df[2] = i_df[D1Idx(i_x, ym)];
            out_df[3] = i_df[D1Idx(xp, i_y)];
            out_df[4] = i_df[D1Idx(i_x, yp)];
            out_df[5] = i_df[D1Idx(xm, ym)];
            out_df[6] = i_df[D1Idx(xp, ym)];
            out_df[7] = i_df[D1Idx(xp, yp)];
            out_df[8] = i_df[D1Idx(xm, yp)];
            return out_df;
        }

        private (double[], double, double[]) UpdateStep(double[,] i_df, int i_x, int i_y)
        {
            Debug.Log("UpdateStep");
            double[] loc_df = Streaming(i_df.Cast<double>().ToArray(), i_x, i_y);
            int mapID = WallMap[D1Idx(i_x, i_y)];
            double loc_density;
            double[] loc_velocity;
            if (mapID == wallID)
            {
                loc_df = Bounceback(loc_df);
                loc_density = 1;
                loc_velocity = new double[] { 0.0, 0.0 };
            }
            else if (mapID == inflowID)
            {
                (loc_df, loc_density, loc_velocity) = InflowBoundary(i_x, i_y, i_df.Cast<double>().ToArray());
            }
            else if (mapID == outflowID)
            {
                loc_df = OutflowCondition(i_x, i_y, loc_df);
                (loc_density, loc_velocity) = ComputeMacro(loc_df);
            }
            else
            {
                (loc_density, loc_velocity) = ComputeMacro(loc_df);
                loc_df = loc_df.Select((df, idx) => df + omega * (Equilibrium(loc_density, loc_velocity[0], loc_velocity[1])[idx] - df)).ToArray();
            }
            return (loc_df, loc_density, loc_velocity);
        }

        public void NextIteration()
        {
            Debug.Log("NextIteration");
            int t = currentIteration;
            if (t < numberOfIterations)
            {
                currentIteration++;
                for (int x = 0; x < nx; x++)
                {
                    for (int y = 0; y < ny; y++)
                    {
                        int idx = D1Idx(x, y);
                        if (t % 2 == 0)
                        {
                            var (new_df, new_density, new_velocity) = UpdateStep(df1, x, y);
                            for (int i = 0; i < 9; i++)
                            {
                                df2[idx, i] = new_df[i];
                            }
                            density[idx] = new_density;
                            velocity[idx, 0] = new_velocity[0];
                            velocity[idx, 1] = new_velocity[1];
                        }
                        else
                        {
                            var (new_df, new_density, new_velocity) = UpdateStep(df2, x, y);
                            for (int i = 0; i < 9; i++)
                            {
                                df1[idx, i] = new_df[i];
                            }
                            density[idx] = new_density;
                            velocity[idx, 0] = new_velocity[0];
                            velocity[idx, 1] = new_velocity[1];
                        }
                    }
                }
            }
        }

        public double PhysicalTime()
        {
            Debug.Log("PhysicalTime");
            return currentIteration * dt;
        }

        public int D1Idx(int i_x, int i_y)
        {
            Debug.Log("D1Idx");
            return i_x + i_y * nx;
        }

        private (double, double) PositionCellCenter(int i_x, int i_y)
        {
            Debug.Log("PositionCellCenter");
            return (i_x * dx, i_y * dx - dx / 2.0);
        }

        private (double, double) PositionCell(int i_x, int i_y)
        {
            Debug.Log("PositionCell");
            return (i_x * dx - dx / 2.0, i_y * dx - dx);
        }

        private double Phys2LBMVelocity(double i_u)
        {
            Debug.Log("Phys2LBMVelocity");
            return i_u * dt / dx;
        }

        public double LBM2PhysVelocity(double i_u)
        {
            Debug.Log("LBM2PhysVelocity");
            return i_u * dx / dt;
        }

        private void MapInit()
        {
            Debug.Log("MapInit");
            for (int x = 0; x < nx; x++)
            {
                int idxTop = D1Idx(x, ny - 1);
                int idxDown = D1Idx(x, 0);
                WallMap[idxTop] = wallID;
                WallMap[idxDown] = wallID;
            }

            for (int y = 1; y < ny - 1; y++)
            {
                int idxLeft = D1Idx(0, y);
                int idxRight = D1Idx(nx - 1, y);
                WallMap[idxLeft] = inflowID;
                WallMap[idxRight] = outflowID;
            }

            for (int x = 0; x < nx; x++)
            {
                for (int y = 0; y < ny; y++)
                {
                    var (xf, yf) = PositionCellCenter(x, y);
                    if (Math.Pow(xf - Cx, 2) + Math.Pow(yf - Cy, 2) <= Math.Pow(R, 2))
                    {
                        int idx1d = D1Idx(x, y);
                        WallMap[idx1d] = wallID;
                    }
                }
            }
        }

        private void MacroInit()
        {
            Debug.Log("MacroInit");
            for (int x = 0; x < nx; x++)
            {
                for (int y = 0; y < ny; y++)
                {
                    int idx1d = D1Idx(x, y);
                    velocity[idx1d, 0] = Phys2LBMVelocity(0.0);
                    velocity[idx1d, 1] = Phys2LBMVelocity(0.0);
                    density[idx1d] = 1.0;
                }
            }
        }

        private void DistributionInit()
        {
            Debug.Log("DistributionInit");
            for (int x = 0; x < nx; x++)
            {
                for (int y = 0; y < ny; y++)
                {
                    int idx1d = D1Idx(x, y);
                    var eq = Equilibrium(density[idx1d], velocity[idx1d, 0], velocity[idx1d, 1]);
                    for (int i = 0; i < 9; i++)
                    {
                        df1[idx1d, i] = eq[i];
                    }
                }
            }
        }

        public void Start()
        {
            Debug.Log("LBM Start");
            MapInit();
            MacroInit();
            DistributionInit();
        }

        private (double, double) InflowVelocity(double i_y)
        {
            Debug.Log("InflowVelocity");
            double u0 = Phys2LBMVelocity(uMax);
            return (4.0 * u0 * i_y * (ny * dx - i_y) / Math.Pow(ny * dx, 2), 0.0);
        }

        private (double, double[]) ComputePhysMacro(double[] i_df)
        {
            Debug.Log("ComputePhysMacro");
            double nonDimCS2 = Phys2LBMVelocity(Phys2LBMVelocity(1.0 / 3.0));
            double rho_loc = i_df.Sum();
            double vx_loc = i_df.Select((df, idx) => df * velocitySet[0, idx]).Sum() / rho_loc;
            double vy_loc = i_df.Select((df, idx) => df * velocitySet[1, idx]).Sum() / rho_loc;
            return ((rho_loc - 1.0) * cs2 + p0, new double[] { LBM2PhysVelocity(vx_loc), LBM2PhysVelocity(vy_loc) });
        }
    }
}