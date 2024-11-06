using UnityEngine;
using UnityEngine.Video;

public class VideoController : MonoBehaviour
{
    private VideoPlayer videoPlayer;
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        videoPlayer = GetComponent<VideoPlayer>();
        videoPlayer.Pause(); // Ensure video starts paused
    }

    public void PlayVideo()
    {
        if (videoPlayer.isPlaying)
        {
            videoPlayer.Pause();
        }
        else
        {
            videoPlayer.Play();
        }
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
