# free_moving_mice
Mouse movement analysis scripts developed for the Biophotonics Imaging Technology (BIT) Lab. Specifically:

- plot speed and angular speed frequency histogram
- plot speed and angular speed over time
- plot trajectory color-coded by speed or angular speed
- create movie of mouse movement

## Usage
The following changes can be made in *free_moving_mice.ipynb*

Read in a file:
> df = pd.read_csv("file_name.csv").iloc[2:]

Change constants:
> LPP = new_value  
TPF = new_value

Save a vector to csv:
> pd.DataFrame(vector).to_csv("vector_name.csv", index=None)

Adjust threshold for trajectory plot:
> speed_t[speed_t > threshold] = threshold  
angular_speed_t[angular_speed_t > threshold] = threshold

## Mouse Movement Movie
https://user-images.githubusercontent.com/98730743/197016815-2e6ab83e-b37a-4bb8-a337-1d2a40d75559.mp4

