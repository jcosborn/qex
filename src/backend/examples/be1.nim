import backend/accel

block:
  var x = 1.0'f32
  let y = 2.0'f32
  const z = 3.0'f32
  echo "x: ", x
  onGpu:
    x = y * z
  echo "x: ", x
