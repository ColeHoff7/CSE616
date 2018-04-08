echo Converting images to GIF
magick.exe convert -delay 0.5 H*.jpg plot.gif
echo Collecting garbage
rm *.jpg*