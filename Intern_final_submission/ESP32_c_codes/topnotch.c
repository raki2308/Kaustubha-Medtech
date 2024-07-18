#include <WiFi.h>
#include <BluetoothSerial.h>

const int analogInputPin = A0;  // Analog input pin
const int outputPin = 9;       // Digital output pin

int sensorValue = 0;           // Variable to store the raw sensor value
float filteredValue = 0.0;     // Variable to store the filtered value

// Filter parameters
const float samplingFrequency = 1000.0;  // Adjust this based on your actual sampling rate
const float centerFrequency = 50.0;      // Desired center frequency in Hz
const float bandwidth = 2.0;             // Bandwidth in Hz

// Coefficients for a 2nd order notch filter
const float Q = centerFrequency / bandwidth;
const float omega = 2.0 * PI * centerFrequency / samplingFrequency;
const float alpha = sin(omega) / (2.0 * Q);
const float a0 = 1.0 + alpha;
const float a1 = -2.0 * cos(omega);
const float a2 = 1.0 - alpha;
const float b0 = 1.0;
const float b1 = -2.0 * cos(omega);
const float b2 = 1.0;

BluetoothSerial SerialBT;       // Bluetooth serial instance

const char* ssid = "your-ssid";      // Replace with your Wi-Fi SSID
const char* password = "your-password";  // Replace with your Wi-Fi password

void setup() {
  pinMode(outputPin, OUTPUT);
  Serial.begin(9600);

  // Initialize Wi-Fi
  WiFi.begin(ssid, password);
  while (WiFi.status() != WL_CONNECTED) {
    delay(1000);
    Serial.println("Connecting to WiFi...");
  }
  Serial.println("Connected to WiFi");

  // Initialize Bluetooth
  SerialBT.begin("ESP32_BT");  // Bluetooth device name

  // Print the IP address of the ESP32
  Serial.print("IP Address: ");
  Serial.println(WiFi.localIP());
}

void loop() {
  // Read the raw sensor value
  sensorValue = analogRead(analogInputPin);

  // Apply the 2nd order notch filter
  filteredValue = notchFilter(sensorValue, filteredValue);

  // Output the filtered value to the digital pin
  analogWrite(outputPin, filteredValue / 4);  // Adjust the division factor as needed

  // Print the raw and filtered values to the serial monitor
  Serial.print("Raw: ");
  Serial.print(sensorValue);
  Serial.print("   Filtered: ");
  Serial.println(filteredValue);

  // Send data over Bluetooth
  SerialBT.print("Raw: ");
  SerialBT.print(sensorValue);
  SerialBT.print("   Filtered: ");
  SerialBT.println(filteredValue);

  delay(10);  // Adjust the delay as needed
}

float notchFilter(float input, float outputPrev) {
  static float x[3] = {0.0, 0.0, 0.0};
  static float y[3] = {0.0, 0.0, 0.0};

  x[0] = x[1];
  x[1] = x[2];
  x[2] = input / 4095.0;  // Normalize input to the range [0, 1] for ESP32

  y[0] = y[1];
  y[1] = y[2];

  y[2] = (b0 * x[2] + b1 * x[1] + b2 * x[0] - a1 * y[1] - a2 * y[0]) / a0;

  return y[2];
}
