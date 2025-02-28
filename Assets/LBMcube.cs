using System;
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using System.Threading;
using UnityEngine;
using LBMNamespace;

public class LBMcube : MonoBehaviour
{
    public Texture2D texture;
    public Vector2 textureSize = new Vector2(2048, 2048);
    private LBM lbm;
    private Thread simulationThread;
    private bool isRunning = true;
    private bool updateTexture = false;

    void Start()
    {
        Debug.Log("LBMcube Start");
        var r = GetComponent<Renderer>();
        texture = new Texture2D((int)textureSize.x, (int)textureSize.y);
        r.material.mainTexture = texture;

        // Initialize LBM simulation
        double uMax = 0.1;
        double nu_phys = 0.01;
        double Lx = 2.0;
        double Ly = 2.0;
        double T = 10.0;
        double dx = 0.01;
        lbm = new LBM(uMax, nu_phys, Lx, Ly, T, dx);
        lbm.Start();

        // Start the simulation with a delay
        StartCoroutine(StartSimulationWithDelay(2.0f));
    }

    void OnDestroy()
    {
        Debug.Log("LBMcube OnDestroy");
        // Ensure the simulation thread is properly terminated
        isRunning = false;
        if (simulationThread != null && simulationThread.IsAlive)
        {
            simulationThread.Join();
        }
    }

    IEnumerator StartSimulationWithDelay(float delay)
    {
        Debug.Log("StartSimulationWithDelay");
        yield return new WaitForSeconds(delay);
        simulationThread = new Thread(RunSimulation);
        simulationThread.Start();
    }

    void RunSimulation()
    {
        Debug.Log("RunSimulation started");
        int totalIterations = 100;

        for (int i = 0; i < totalIterations && isRunning; i++)
        {
            lbm.NextIteration();
            updateTexture = true;
            Debug.Log("Iteration " + i + " completed, updateTexture set to true");
            Thread.Sleep(10); // Sleep to simulate work and allow other threads to run
        }
        Debug.Log("RunSimulation completed");
    }

    void Update()
    {
        if (updateTexture)
        {
            Debug.Log("Update called, updating texture");
            UpdateTexture();
            updateTexture = false;
        }
    }

    void UpdateTexture()
{
    Debug.Log("Updating texture...");
    Color[] pixels = new Color[texture.width * texture.height];

    for (int x = 0; x < texture.width; x++)
    {
        for (int y = 0; y < texture.height; y++)
        {
            int idx = lbm.D1Idx(x, y);
            double[] df = new double[9];
            for (int i = 0; i < 9; i++)
            {
                df[i] = lbm.df1[idx, i];
            }
            var (rho_loc, velocity) = lbm.ComputeMacro(df);
            double vx = velocity[0];
            double vy = velocity[1];

            // Map density and velocity to heatmap color
            float magnitude = Mathf.Clamp((float)Math.Sqrt(vx * vx + vy * vy), 0, 1);
            Color color = Color.Lerp(Color.blue, Color.red, magnitude);
            pixels[x + y * texture.width] = color;
        }
    }

    texture.SetPixels(pixels);
    texture.Apply();
    Debug.Log("Texture updated.");
}
}