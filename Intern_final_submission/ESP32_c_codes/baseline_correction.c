#include <WiFi.h>
#include <BluetoothSerial.h>

const int analogInputPin = A0;  // Analog input pin (replace A0 with the appropriate GPIO number)
const int outputPin = 9;       // Digital output pin (replace 9 with the appropriate GPIO number)

int sensorValue = 0;           // Variable to store the raw sensor value
float baseline = 0.0;          // Variable to store the baseline value
float correctedValue = 0.0;    // Variable to store the baseline-corrected value

// Filter parameters for moving average baseline correction
const int numSamples = 10;      // Number of samples to average
int samples[numSamples];        // Array to store past samples
int sampleIndex = 0;            // Index for the sample array
float sum = 0.0;                // Running sum of past samples

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

  // Update the baseline using a moving average
  updateBaseline(sensorValue);

  // Perform baseline correction
  correctedValue = sensorValue - baseline;

  // Output the corrected value to the digital pin
  analogWrite(outputPin, correctedValue / 4);  // Adjust the division factor as needed

  // Print the raw, baseline, and corrected values to the serial monitor
  Serial.print("Raw: ");
  Serial.print(sensorValue);
  Serial.print("   Baseline: ");
  Serial.print(baseline);
  Serial.print("   Corrected: ");
  Serial.println(correctedValue);

  // Send data over Bluetooth
  SerialBT.print("Raw: ");
  SerialBT.print(sensorValue);
  SerialBT.print("   Baseline: ");
  SerialBT.print(baseline);
  SerialBT.print("   Corrected: ");
  SerialBT.println(correctedValue);

  delay(10);  // Adjust the delay as needed
}
 
void updateBaseline(int newValue) {
  // Subtract the oldest sample from the running sum
  sum -= samples[sampleIndex];

  // Add the new sample to the sum
  sum += newValue;

  // Store the new sample in the array
  samples[sampleIndex] = newValue;

  // Increment the sample index and wrap around if necessary
  sampleIndex = (sampleIndex + 1) % numSamples;

  // Calculate the new baseline value
  baseline = sum / numSamples;
}
