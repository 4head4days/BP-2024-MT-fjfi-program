using System;
using System.IO;
using System.Net.Sockets;
using System.Runtime.InteropServices;
using UnityEngine;

public class SocketImageReceiver : MonoBehaviour
{
    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    public struct Pixel
    {
        public byte r;
        public byte g;
        public byte b;
    }

    private const string ServerIP = "127.0.0.1";
    private const int ServerPort = 12345;

    private Texture2D texture;

    void Start()
    {
        try
        {
            // Connect to the simulation server
            using (TcpClient client = new TcpClient(ServerIP, ServerPort))
            using (NetworkStream stream = client.GetStream())
            {
                Debug.Log("Connected to simulation application");

                // Read the dimensions (8 bytes: 2 integers)
                byte[] dimensionsBuffer = new byte[8];
                ReadFully(stream, dimensionsBuffer, dimensionsBuffer.Length);

                int width = BitConverter.ToInt32(dimensionsBuffer, 0);
                int height = BitConverter.ToInt32(dimensionsBuffer, 4);

                Debug.Log($"Received dimensions: {width}x{height}");

                // Create the texture
                texture = new Texture2D(width, height, TextureFormat.RGB24, false);

                // Read the pixel data (width * height * 3 bytes)
                int pixelCount = width * height;
                byte[] pixelBuffer = new byte[pixelCount * 3];
                ReadFully(stream, pixelBuffer, pixelBuffer.Length);

                Debug.Log($"Received pixel data: {pixelBuffer.Length} bytes");

                // Populate the texture
                int bufferIndex = 0;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        byte r = pixelBuffer[bufferIndex++];
                        byte g = pixelBuffer[bufferIndex++];
                        byte b = pixelBuffer[bufferIndex++];
                        texture.SetPixel(x, y, new Color(r / 255f, g / 255f, b / 255f));
                    }
                }
                texture.Apply();

                // Assign the texture to a plane
                Renderer renderer = GetComponent<Renderer>();
                if (renderer != null)
                {
                    renderer.material.mainTexture = texture;
                }
                else
                {
                    Debug.LogWarning("No Renderer found on the GameObject");
                }

                Debug.Log("Texture updated with simulation data");
            }
        }
        catch (SocketException ex)
        {
            Debug.LogError($"SocketException: {ex.Message}");
        }
        catch (IOException ex)
        {
            Debug.LogError($"IOException: {ex.Message}");
        }
        catch (Exception ex)
        {
            Debug.LogError($"Unexpected error: {ex.Message}");
        }
    }
    
    private void ReadFully(NetworkStream stream, byte[] buffer, int size)
    {
        int bytesRead = 0;
        while (bytesRead < size)
        {
            int read = stream.Read(buffer, bytesRead, size - bytesRead);
            if (read == 0)
            {
                throw new IOException("Connection closed while reading data");
            }
            bytesRead += read;
        }
    }
}
