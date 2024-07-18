#include <WiFi.h>
#include <BluetoothSerial.h>

const int analogInputPin = A0;  // Analog input pin (replace A0 with the appropriate GPIO number)
const int outputPin = 9;       // Digital output pin (replace 9 with the appropriate GPIO number)

int sensorValue = 0;           // Variable to store the raw sensor value
float filteredValue = 0.0;     // Variable to store the filtered value

// Filter parameters
const float samplingFrequency = 1000.0;  // Adjust this based on your actual sampling rate
const float cutoffFrequency = 0.33;      // Desired cutoff frequency in Hz

// Coefficients for a 2nd order high-pass Butterworth filter
const float b0 = 0.963174;  // 1 - cos(2 * PI * cutoff / samplingFrequency)
const float b1 = -1.926349; // 2 * cos(2 * PI * cutoff / samplingFrequency)
const float b2 = 0.963174;  // 1 - cos(2 * PI * cutoff / samplingFrequency)
const float a1 = -1.911402; // -2 * cos(2 * PI * cutoff / samplingFrequency)
const float a2 = 0.914975;  // cos(2 * PI * cutoff / samplingFrequency) - 1

BluetoothSerial SerialBT;  // Bluetooth serial instance

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

  // Apply the 2nd order high-pass Butterworth filter
  filteredValue = highPassButterworthFilter(sensorValue, filteredValue);

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

float highPassButterworthFilter(float input, float outputPrev) {
  static float x[3] = {0.0, 0.0, 0.0};
  static float y[3] = {0.0, 0.0, 0.0};

  x[0] = x[1];
  x[1] = x[2];
  x[2] = input / 4095.0;  // Normalize input to the range [0, 1]

  y[0] = y[1];
  y[1] = y[2];

  y[2] = b0 * x[2] + b1 * x[1] + b2 * x[0] + a1 * y[1] + a2 * y[0];

  return y[2];
}
