using UnityEngine;

public class ExpandingCircleAnimation : MonoBehaviour
{
    public float animationSpeed = 1.0f; // Speed of the animation
    public int textureResolution = 256; // Size of the texture (e.g., 256x256)
    private Texture2D animatedTexture;
    private Material material;

    private float radius; // Current radius of the circle

    void Start()
    {
        // Create a new texture
        animatedTexture = new Texture2D(textureResolution, textureResolution);
        animatedTexture.wrapMode = TextureWrapMode.Clamp;

        // Get the material of the box and assign the texture
        material = GetComponent<Renderer>().material;
        material.mainTexture = animatedTexture;
    }

    void Update()
    {
        // Increase the radius over time, and reset it when it exceeds the texture's bounds
        radius += Time.deltaTime * animationSpeed * (textureResolution / 2); // Scale with texture size
        if (radius > textureResolution / 2)
        {
            radius = 0;
        }

        // Generate the expanding circle texture
        GenerateExpandingCircleTexture();

        // Apply changes to the texture
        animatedTexture.Apply();
    }

    void GenerateExpandingCircleTexture()
    {
        // Clear the texture (set all pixels to black)
        for (int x = 0; x < textureResolution; x++)
        {
            for (int y = 0; y < textureResolution; y++)
            {
                animatedTexture.SetPixel(x, y, Color.black);
            }
        }

        // Draw the circle
        for (int x = 0; x < textureResolution; x++)
        {
            for (int y = 0; y < textureResolution; y++)
            {
                // Calculate the distance from the center
                float distance = Vector2.Distance(new Vector2(x, y), new Vector2(textureResolution / 2, textureResolution / 2));

                // Draw a white pixel if the distance is close to the current radius
                if (Mathf.Abs(distance - radius) < 1.5f) // Controls circle thickness
                {
                    animatedTexture.SetPixel(x, y, Color.white);
                }
            }
        }
    }
}
