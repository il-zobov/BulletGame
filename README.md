# BulletGame
Test app for 4A games studio - small 2D shooting game engine example.

A list of segments (walls), represented by two points on a plane, is predefined.

Your task is to write BulletManager class, consisting of two public methods:
- void Update (float time), where time – global update time in seconds. This method calculates
bullet movement in given time, and in case of collision with the wall, removes the wall from the
list and bullet continues movement with its trajectory reflected.

- void Fire (float2 pos, float2 dir, float speed, float time, float life_time), where pos – starting
point of the bullet on the plane (in meters), dir – direction, speed – bullet speed (m/s), time –
when the bullet will be fired, life_time – time before bullet disappears. This method adds bullet
to manager for further updates. Update method must be called from main thread, but Fire
method could be called from different threads.
