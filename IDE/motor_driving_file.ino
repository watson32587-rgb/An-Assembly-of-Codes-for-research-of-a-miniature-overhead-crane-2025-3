// โค้ดทดสอบ Hardware เบื้องต้น
const int stepPin = 9;
const int dirPin = 8;

void setup() {
  pinMode(stepPin, OUTPUT);
  pinMode(dirPin, OUTPUT);
  digitalWrite(dirPin, LOW); // กำหนดทิศทาง
}

void loop() {
  // สั่งให้หมุนด้วยความเร็วคงที่
  digitalWrite(stepPin, HIGH);
  delayMicroseconds(469); // ปรับเลขนี้: น้อย = เร็ว, มาก = ช้า
  digitalWrite(stepPin, LOW);
  delayMicroseconds(469);
}