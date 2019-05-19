# PhysicsAABB
This project demonstrates some physics done without bullet, and uses AABBs

By: Angus Poole

## About
Build/Compile in Release/x86.

This project uses nothing but math to perform physics - no bullet involved. Acceleration forces are applied to the ship to move, and because the ship is in space, it continues to drift afterwards.

Rotating the ship changes the direction new forces are applied, but not existing forces.

The terrain is divided using AABBs, and whatever boxes the ship is inside are visible.

There are 3 different modes: basic, collision spheres, and AABB lines. The basic mode is the default one. The collision sphere mode shows the spheres used for collisions on the ship. The AABB lines mode shows lines from the min xyz to the max xyz of each box in the collection of AABBs.

If the ship collides with the terrain, a smoke particle effect should spawn on the terrain at the collision point.


## Controls
* Use WASD to apply thrust forces to the ship
* Use IJKL to rotate the ship
* Press 1 to go to "basic" mode
* Press 2 to go to "collision spheres" mode
* Press 3 to go to "AABB lines" mode