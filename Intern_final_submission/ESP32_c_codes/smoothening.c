#include <WiFi.h>
#include <BluetoothSerial.h>

const int analogInputPin = A0;  // Analog input pin
const int outputPin = 9;       // Digital output pin

int ecgValue = 0;              // Variable to store the raw ECG signal
float smoothedValue = 0.0;     // Variable to store the smoothed ECG signal

// Filter parameters for moving average smoothing
const int numSamples = 10;      // Number of samples to average
float samples[numSamples];      // Array to store past samples
int sampleIndex = 0;            // Index for the sample array

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
  // Read the raw ECG signal
  ecgValue = analogRead(analogInputPin);

  // Smooth the ECG signal using a moving average
  smoothedValue = movingAverageFilter(ecgValue);

  // Output the smoothed ECG value to the digital pin
  analogWrite(outputPin, smoothedValue / 16);  // Adjust the division factor as needed for ESP32

  // Print the raw and smoothed ECG values to the serial monitor
  Serial.print("Raw ECG: ");
  Serial.print(ecgValue);
  Serial.print("   Smoothed ECG: ");
  Serial.println(smoothedValue);

  // Send data over Bluetooth
  SerialBT.print("Raw ECG: ");
  SerialBT.print(ecgValue);
  SerialBT.print("   Smoothed ECG: ");
  SerialBT.println(smoothedValue);

  delay(10);  // Adjust the delay as needed
}

float movingAverageFilter(int newValue) {
  // Add the new sample to the array
  samples[sampleIndex] = newValue;

  // Increment the sample index and wrap around if necessary
  sampleIndex = (sampleIndex + 1) % numSamples;

  // Calculate the moving average
  float sum = 0.0;
  for (int i = 0; i < numSamples; ++i) {
    sum += samples[i];
  }

  return sum / numSamples;
}
