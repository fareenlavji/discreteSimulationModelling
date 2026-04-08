# 1. Overview
This project involves simulating a computer processing system composed of two processing nodes, multiple buffers, and a router. The goal is to analyze system performance using discrete-event simulation and statistical analysis with 95% confidence intervals.

---

# 2. System Description

## 2.1. Arrival Processes

### 2.1.1. Buffer 1
- **Packet type**: Type I
- **Arrival process**: Poisson
- **Mean interarrival time**: 4 ms

### 2.1.2. Buffer 2
- **Packet type**: Type II
- **Arrival process**: Poisson
- **Mean interarrival time**: 12 ms

## 2.2. Processing Node 1
### 2.2.1. Service Policy
Alternates between:
- Buffer 1 for 50 ms
- Buffer 2 for 30 ms

If one buffer is empty, the processor serves the non-empty buffer continuously. Once service of a packet starts, it must complete before switching buffers.

### 2.2.2. Service Time Distributions
- Buffer 1 (Type I): Uniform (1 ms, 3 ms)
- Buffer 2 (Type II): Uniform (2 ms, 6 ms)

## 2.3. Routing After Processing Node 1

### 2.3.1. Type I packets
Routed to Processing Node 2 with probability 1

### 2.3.2. Type II packets
- Routed to Processing Node 2 with probability 0.5
- Otherwise discarded from the system

## 2.4. Router 1 Operation
- Redirects a packet away from the path if more than 5 packets are waiting ahead in the buffer
- Routing operation time is assumed to be zero
- Used to compute packet redirection probability

## 2.5. Processing Node 2
- Two parallel processors
- Shared common buffer
- Processing time at each processor:
    - Exponentially distributed
    - Mean = 5 ms

---

# 3. Quantities of Interest
We are required to estimate the following (each with 95% confidence intervals, non-simultaneous):
- [ ] Average number of Type I packets waiting and being served in Buffer 1 / Processor 1
- [ ] Average number of Type II packets waiting and being served in Buffer 2 / Processor 1
- [ ] Average number of Type I packets that must wait in Buffer 1
- [ ] Average number of Type II packets that must wait in Buffer 2
- Average waiting time in:
    - [ ] Buffer 1
    - [ ] Buffer 2
    - [ ] Include histograms
- Average total time (waiting + service) at Processor 1 for:
    - [ ] Type I packets
    - [ ] Type II packets
    - [ ] Include histograms
- [ ] Distribution fitting: Test whether results from (5) and (6) match known theoretical distributions
- Processor 2 analysis:
    - [ ] Histogram of waiting time (buffer + processor)
    - [ ] Test for known distribution fit
- [ ] Probability that a packet is redirected by Router 1

---

# 4. Implementation Plan
We need to integrate the 12-step simulation project cycle with the specific technical requirements of the system.

## 4.1. Phase 1: Model Conceptualization & Logic Design
Outline of how the simulation will track time and state changes.

### 4.1.1. Defined State Variables
- `Clock`: The current simulation time.
- `Buffer_1_Count`, `Buffer_2_Count`: Number of packets in Node 1 buffers.
- `Node_2_Buffer_Count`: Number of packets waiting for the two processors at Node 2.
- `Node_1_Status`: (Idle/Busy with Buffer 1/Busy with Buffer 2).
- `Processor_Status`: (Idle/Busy) for each of the two processors at Node 2.

### 4.1.2. Identified Primary Events
- **Arrival (Type I & II)**: Schedules the next arrival and updates Buffer counts.
- **Departure from Node 1**: Triggers routing logic (discard, redirect, or send to Node 2) and schedules the next Node 1 departure if the current buffer is still being served.
- **Departure from Node 2**: Frees a processor and pulls the next packet from the common buffer.
- **Node 1 Buffer Switch**: A "timer" event to handle the alternation logic (50ms for Buffer 1, 30ms for Buffer 2).

## 4.2. Phase 2: Input Modeling (Random Variate Generation)
Implement mathematical algorithms to generate stochastic values for known distributions.

### 4.2.1. Interarrival Times (Inverse-Transform)
- Type I: $X = -4 \ln(U)$.
- Type II: $X = -12 \ln(U)$.

### 4.2.2. Service Times (Uniform/ Exponential)
- Node 1 (Buffer 1): $X = 1 + (3 - 1)U$.
- Node 1 (Buffer 2): $X = 2 + (6 - 2)U$.
- Node 2 Processors: $X = -5 \ln(U)$.

### 4.2.3. Stochastic Decisions
Use a random number $U \sim U(0,1)$ to decide if a Type II packet is discarded ($U \le 0.5$).

## 4.3. Phase 3: Implementation & Programming
Choose a high-level language (Python, Java, or MATLAB) and build the discrete-event engine.

### 4.3.1. The Future Event List (FEL)
Use a priority queue to store and sort events by their scheduled time.

### 4.3.2. The Simulation Loop
1. Fetch the earliest event from the FEL and advance the `Clock`.
2. Execute the corresponding event subroutine (Arrival, Departure, or Switch).
3. Update **statistical accumulators** (sum of wait times, total packets redirected, etc.).

### 4.3.3. Node 1 Logic Implementation
Ensure code handles the rule that if one buffer is empty, the processor serves the other continuously, and that service must be completed before switching buffers.

### 4.3.4. Router Logic
Before moving a packet to Node 2, check if `Node_2_Buffer_Count > 5`. If so, increment the "redirected" counter.

## 4.4. Phase 4: Production Runs & Output Analysis
After verifying the code is bug-free, perform production runs to generate the required data.
1. **Execution:** Conduct multiple independent replications to ensure statistical significance.
2. **Calculate Point Estimates**:
    a. **Average Waiting Time ($\hat{W}$)**: $\frac{\text{Total Waiting Time}}{\text{Total Packets}}$.
    b. **Average Queue Length ($\hat{L}$)**: $\frac{\sum (\text{Queue Length} \times \text{Duration})}{T}$.
3. **Redirection Probability**: Calculate the ratio of redirected packets to total packets exiting Node 1.
4. **Confidence Intervals**: Calculate a **95% Confidence Interval** for each metric using the $t$-distribution critical value (e.g., $1.96 \times \frac{S}{\sqrt{n}}$ for large samples).

## 4.5. Phase 5: Statistical Matching & Reporting
1. **Histogram Construction**: Create frequency distributions for the waiting times in Buffers 1 and 2, and the waiting time at Node 2.
2. **Goodness-of-Fit**: Use the **Chi-Square test** to see if your output histograms match known distributions (like Exponential).
3. **Final Report**: Document your logic, source code, and statistical results according to the marking scheme provided in your project guidelines.
