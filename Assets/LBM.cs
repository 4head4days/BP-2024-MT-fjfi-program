using System;

namespace LBMNamespace{
class LBM
{
    // Fixed Parameters
    private const double p0 = 0;
    private const double nu_LB = 0.01;
    private const double cs2 = 1.0 / 3.0;
    private const double tau_LB = (nu_LB + 0.5) / cs2;

    // Domain Parameters
    private double Lx = 2.2;
    private double Ly = 0.41;
    private double T = 1000;
    private double nu = 0.001;
    private double u0 = 0.3;
    private double rho_0 = 1;

    // Cylinder Parameters
    private double Cx = 0.2;
    private double Cy = 0.2;
    private double R = 0.05;

    // Grid Resolution
    private double dx = 0.01;
    private int nx;
    private int ny;
    private double dt;
    private int numberOfIterations;
    private int currentIteration;

    private double[] density;
    private double[,] velocity;
    private int[] WallMap;
    private double[,] df1;
    private double[,] df2;

    // Boundary Mapping
    private const int flowID = 0;
    private const int wallID = 1;
    private const int inflowID = 2;
    private const int outflowID = 3;

    // Lattice Parameters
    private double omega;
    private double[] weights = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
    private double[,] velocitySet =
    {
        { 0, 1, 0, -1, 0, 1, -1, -1, 1 },
        { 0, 0, 1, 0, -1, 1, 1, -1, -1 }
    };

    public LBM()
    {
        nx = (int)Math.Ceiling(Lx / dx);
        ny = (int)Math.Ceiling(Ly / dx);
        dt = nu_LB / nu * dx * dx;
        numberOfIterations = (int)Math.Ceiling(T / dt);
        currentIteration = 0;

        density = new double[nx * ny];
        velocity = new double[nx * ny, 2];
        WallMap = new int[nx * ny];
        df1 = new double[nx * ny, 9];
        df2 = new double[nx * ny, 9];
        omega = cs2 / (nu_LB + 0.5);
    }

    private double[] Equilibrium(double i_rho, double i_u, double i_v)
    {
        double[] eq = new double[9];
        for (int idx = 0; idx < 9; idx++)
        {
            double eu = velocitySet[0, idx] * i_u + velocitySet[1, idx] * i_v;
            eq[idx] = i_rho * weights[idx] * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * (i_u * i_u + i_v * i_v));
        }
        return eq;
    }

    private (double, double[]) ComputeMacro(double[] i_df)
    {
        double rho = 0;
        double vx = 0, vy = 0;
        for (int i = 0; i < 9; i++)
        {
            rho += i_df[i];
            vx += i_df[i] * velocitySet[0, i];
            vy += i_df[i] * velocitySet[1, i];
        }
        return (rho, new double[] { vx / rho, vy / rho });
    }

    private int D1Idx(int i_x, int i_y) => i_x + i_y * nx;

    private void MapInit()
    {
        for (int x = 0; x < nx; x++)
        {
            WallMap[D1Idx(x, ny - 1)] = wallID;
            WallMap[D1Idx(x, 0)] = wallID;
        }
        for (int y = 1; y < ny - 1; y++)
        {
            WallMap[D1Idx(0, y)] = inflowID;
            WallMap[D1Idx(nx - 1, y)] = outflowID;
        }
        for (int x = 0; x < nx; x++)
        {
            for (int y = 0; y < ny; y++)
            {
                double xf = x * dx;
                double yf = y * dx;
                if ((xf - Cx) * (xf - Cx) + (yf - Cy) * (yf - Cy) <= R * R)
                {
                    WallMap[D1Idx(x, y)] = wallID;
                }
            }
        }
    }

    private void MacroInit()
    {
        for (int i = 0; i < nx * ny; i++)
        {
            velocity[i, 0] = 0;
            velocity[i, 1] = 0;
            density[i] = 1;
        }
    }

    private void DistributionInit()
    {
        for (int x = 0; x < nx; x++)
        {
            for (int y = 0; y < ny; y++)
            {
                int idx = D1Idx(x, y);
                df1[idx, ..] = Equilibrium(density[idx], velocity[idx, 0], velocity[idx, 1]);
            }
        }
    }

    public void Start()
    {
        MapInit();
        MacroInit();
        DistributionInit();
    }

    public void NextIteration()
    {
        if (currentIteration < numberOfIterations)
        {
            currentIteration++;
            for (int x = 0; x < nx; x++)
            {
                for (int y = 0; y < ny; y++)
                {
                    int idx = D1Idx(x, y);
                    if (currentIteration % 2 == 0)
                    {
                        df2[idx, ..] = df1[idx, ..];
                    }
                    else
                    {
                        df1[idx, ..] = df2[idx, ..];
                    }
                }
            }
        }
    }
}
}