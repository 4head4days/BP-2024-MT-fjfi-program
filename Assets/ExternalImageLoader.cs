using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class ExternalImageLoader : MonoBehaviour
{
    public string imageFolderPath = @"D:\Blbosti reborn\Blbosti"; // Set the folder path here
    public float displayTime = 1.0f; // Time per image in seconds
    private Material material;
    private List<Texture2D> loadedImages = new List<Texture2D>();
    private int currentImageIndex = 0;

    void Start()
    {
        material = GetComponent<Renderer>().material;

        // Load images from the folder
        if (Directory.Exists(imageFolderPath))
        {
            LoadImagesFromFolder();
            if (loadedImages.Count > 0)
            {
                StartCoroutine(PlayImageSequence());
            }
            else
            {
                Debug.LogError("No images found in the specified folder.");
            }
        }
        else
        {
            Debug.LogError("The specified folder path does not exist: " + imageFolderPath);
        }
    }

    void LoadImagesFromFolder()
    {
        string[] imageFiles = Directory.GetFiles(imageFolderPath, "*.png"); // Load PNG images (add "*.jpg" or others if needed)
        
        foreach (string imageFile in imageFiles)
        {
            Texture2D texture = LoadTextureFromFile(imageFile);
            if (texture != null)
            {
                loadedImages.Add(texture);
            }
        }
    }

    Texture2D LoadTextureFromFile(string filePath)
    {
        byte[] fileData = File.ReadAllBytes(filePath);
        Texture2D texture = new Texture2D(2, 2);
        if (texture.LoadImage(fileData)) // LoadImage will auto-resize the texture
        {
            return texture;
        }
        else
        {
            Debug.LogError("Failed to load texture from: " + filePath);
            return null;
        }
    }

    IEnumerator PlayImageSequence()
    {
        while (true)
        {
            material.mainTexture = loadedImages[currentImageIndex];
            yield return new WaitForSeconds(displayTime);
            currentImageIndex = (currentImageIndex + 1) % loadedImages.Count;
        }
    }
}
