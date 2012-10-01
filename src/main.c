 /**


 *****************************************************************************

 **

 **  File        : main.c

 **

 **  Abstract    : main function.

 **

 **  Functions   : main

 **

 **  Environment : Atollic TrueSTUDIO/STM32

 **                STMicroelectronics STM32F10x Standard Peripherals
 Library

 **

 **  Distribution: The file is distributed “as is,” without any warranty

 **                of any kind.

 **

 **  (c)Copyright Atollic AB.

 **  You may use this file as-is or modify it according to the needs of
 your

 **  project. Distribution of this file (unmodified or modified) is not

 **  permitted. Atollic AB permit registered Atollic TrueSTUDIO(R) users
 the

 **  rights to distribute the assembled, compiled & linked contents of this

 **  file as part of an application binary file, provided that it is built

 **  using the Atollic TrueSTUDIO(R) toolchain.

 **

 **


 *****************************************************************************

 */



 /* Includes */

 #include <stddef.h> /*Internal compiler library */

 #include "stm32F10x.h" /* Includes in the folder/path */

 #include "STM32vldiscovery.h"

 #include <math.h> /* C Math Library */



 /* Private typedef
 -----------------------------------------------------------*/

 #define ADC1_DR_Address    ((u32)0x4001244C)

 /* Private define
 ------------------------------------------------------------*/

 #define DACBUFFERSIZE 250

 #define WAVEFREQ 220 /* 220Hz -- A3 */

 #define WAVEFREQ2 392 /* 392Hz -- G4 */

 #define WAVEFREQ3           450 /* 450Hz  */

 #define TIMER6_PRESCALER 24 /* produces a 1MHz tick */

 #define SYSTEM_CLOCK 24E6 /*STM32 runs at 24MHz */



 /* Private macro
 -------------------------------------------------------------*/



 /* Private variables
 ---------------------------------------------------------*/

 uint16_t DACBuffer[DACBUFFERSIZE]; /* Array for the waveform */

 uint8_t flag = 0; /* Global flag for toggling - use later */

 uint32_t fTimer;

 uint32_t aRRTimer;

 uint32_t otherTIM;

 uint32_t timerFreq;

 uint16_t timerPeriod;

 uint16_t otherPeriod;

 uint16_t G4Period;

 uint16_t n;



 /* Private function prototypes
 -----------------------------------------------*/

 void RCC_Configuration(void);

 void DMA_Configuration(void);

 void NVIC_Configuration(void);

 void GPIO_Configuration(void);

 void UART_Configuration(void);

 void Timer_Configuration(uint16_t wavPeriod, uint16_t preScaler);

 void DAC_Configuration(void);

 void ADC_Configuration(void);


 /*******************************************************************************

 * Function Name  : TIM2_IRQHandler

 * Description    : This function handles TIM2 global interrupt request.

 * Input          : None

 * Output         : None

 * Return         : None


 *******************************************************************************/

 void TIM2_IRQHandler(void)

 {

 //int check = 1;



   if (TIM_GetITStatus(TIM2, TIM_IT_Update) != RESET)

   {

       TIM_ClearITPendingBit(TIM2, TIM_IT_Update);




   }



 }





 /**

   * @brief  Main program.

   * @param  None

   * @retval None

   */



 int main(void)

 {

   /* Calculate the gradient of the Sawtooth */

 // m = (uint16_t) ( 4095 / DACBUFFERSIZE);

 int i = 0;

 uint16_t array[DACBUFFERSIZE];
   /* Calculate frequency of timer */

   fTimer = WAVEFREQ * DACBUFFERSIZE;



   /* G4 timer */

   aRRTimer = WAVEFREQ2 * DACBUFFERSIZE;



   /* other timer */

   otherTIM = WAVEFREQ3 * DACBUFFERSIZE;





 /* Calculate Tick Rate */

 timerFreq = SYSTEM_CLOCK / TIMER6_PRESCALER; /* Timer tick is in Hz */



   /* Calculate period of Timer */

 timerPeriod = (uint16_t)( timerFreq / fTimer );



 /* G4 period */

 G4Period = (uint16_t)(timerFreq/aRRTimer);



 /* OTHER */

 otherPeriod = (uint16_t)(timerFreq / otherTIM);





 /* System Clocks Configuration */

 RCC_Configuration();



 /* NVIC configuration */

 NVIC_Configuration();



 /* Configure the GPIO ports */

 GPIO_Configuration();



 /* Timer Configuration */

 Timer_Configuration( (timerPeriod), TIMER6_PRESCALER );



 /* DAC Configuration */

 ADC_Configuration();

 DAC_Configuration();

 /* DMA Config */
 DMA_Configuration ();




 /* Start ADC1 Software Conversion */
 ADC_SoftwareStartConvCmd(ADC1, ENABLE);


   while (1)

   {
  if(DMA_GetFlagStatus(DMA1_FLAG_TC1))

       for(i = 0;i<DACBUFFERSIZE;i++)
       {
         array[i] = DACBuffer[i];
       }
  /* Do nothing */

  //TIM2_IRQHandler();



   }

 }




 void ADC_Configuration(void)
 {
 ADC_InitTypeDef ADC_InitStructure;

 /* ADC1 configuration
 ------------------------------------------------------*/
 ADC_InitStructure.ADC_Mode = ADC_Mode_Independent;
 ADC_InitStructure.ADC_ScanConvMode = ENABLE;
 ADC_InitStructure.ADC_ContinuousConvMode = ENABLE;
 ADC_InitStructure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_None;
 ADC_InitStructure.ADC_DataAlign = ADC_DataAlign_Right;
 ADC_InitStructure.ADC_NbrOfChannel = 1;
 ADC_Init(ADC1, &ADC_InitStructure);
 /* ADC1 regular channel14 configuration */
 ADC_RegularChannelConfig(ADC1, ADC_Channel_14, 1,
 ADC_SampleTime_239Cycles5);
 //ADC_TempSensorVrefintCmd(ENABLE);


 /* Enable ADC1 external trigger conversion */
 ADC_ExternalTrigConvCmd(ADC1, ENABLE);

 /* Enable ADC1 DMA */
 ADC_DMACmd(ADC1, ENABLE);
 /* Enable ADC1 */
 ADC_Cmd(ADC1, ENABLE);
 /* Enable ADC1 reset calibaration register */
 ADC_ResetCalibration(ADC1);
 /* Check the end of ADC1 reset calibration register */
 while(ADC_GetResetCalibrationStatus(ADC1));
 /* Start ADC1 calibaration */
 ADC_StartCalibration(ADC1);
 /* Check the end of ADC1 calibration */
 while(ADC_GetCalibrationStatus(ADC1));
 }
 /**

   * @brief  Configures the different system clocks.

   * @param  None

   * @retval : None

   */

 void RCC_Configuration(void)

 {

 RCC_ADCCLKConfig(RCC_PCLK2_Div2);

   /* DMA1 clock enable */

   RCC_AHBPeriphClockCmd(RCC_AHBPeriph_DMA1, ENABLE);



   /* TIM2 clock, TIM6 and DAC */

   RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM2 | RCC_APB1Periph_TIM6 |
 RCC_APB1Periph_DAC , ENABLE);


   /* adc*/

     RCC_APB2PeriphClockCmd(RCC_APB2Periph_ADC1, ENABLE);


   /* Enable GPIO A and AFIO (Alternate Function) clocks */

   RCC_APB2PeriphClockCmd( RCC_APB2Periph_GPIOA  | RCC_APB2Periph_AFIO,
 ENABLE);


   /* Enable GPIO A and AFIO (Alternate Function) clocks */

    RCC_APB2PeriphClockCmd( RCC_APB2Periph_GPIOC  | RCC_APB2Periph_AFIO,
 ENABLE);


 }



 /**

   * @brief  Configures the nested vectored interrupt controller.

   * @param  None

   * @retval : None

   */

 void NVIC_Configuration(void)

 {

 NVIC_InitTypeDef NVIC_InitStructure;



 /* Enable the TIM2 global Interrupt */

 NVIC_InitStructure.NVIC_IRQChannel = TIM2_IRQn;

 NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = 0;

 NVIC_InitStructure.NVIC_IRQChannelSubPriority = 1;

 NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;

 NVIC_Init(&NVIC_InitStructure);

 }





 /**

   * @brief  Configures the different GPIO ports.

   * @param  None

   * @retval : None

   */

 void GPIO_Configuration(void)

 {

 GPIO_InitTypeDef GPIO_InitStructure;



 /* Configure DAC GPIO as alternate function push-pull */

 GPIO_InitStructure.GPIO_Pin = GPIO_Pin_4;

 GPIO_InitStructure.GPIO_Speed = GPIO_Speed_50MHz;

 GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AF_PP;

 GPIO_Init(GPIOA, &GPIO_InitStructure);

 /* Configure DAC GPIO as alternate function push-pull */

 GPIO_InitStructure.GPIO_Pin = GPIO_Pin_4;

 //GPIO_InitStructure.GPIO_Speed = GPIO_Speed_50MHz;

 GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AIN;

 GPIO_Init(GPIOC, &GPIO_InitStructure);


 }





 /**

   * @brief  Configures the DMA.

   * @param  None

   * @retval : None

   */

 void DMA_Configuration(void)

 {

   DMA_InitTypeDef DMA_InitStructure;



   /* DMA1 Channel3 (triggered by DAC) Config */

   DMA_DeInit(DMA1_Channel1);

   DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)ADC1_DR_Address;

   DMA_InitStructure.DMA_MemoryBaseAddr = (uint32_t)&DACBuffer;

   DMA_InitStructure.DMA_DIR = DMA_DIR_PeripheralSRC;

   DMA_InitStructure.DMA_BufferSize = DACBUFFERSIZE;

   DMA_InitStructure.DMA_PeripheralInc = DMA_PeripheralInc_Disable;

   DMA_InitStructure.DMA_MemoryInc = DMA_MemoryInc_Enable;

   DMA_InitStructure.DMA_PeripheralDataSize =
 DMA_PeripheralDataSize_HalfWord;/* We are transferring 12 bit values to the
 DAC */

   DMA_InitStructure.DMA_MemoryDataSize = DMA_MemoryDataSize_HalfWord; /*
 We are transferring 12 bit values to the DAC */

   DMA_InitStructure.DMA_Mode = DMA_Mode_Circular; /* Circular means the
 Buffer will wrap on itself */

   DMA_InitStructure.DMA_Priority = DMA_Priority_VeryHigh;

   DMA_InitStructure.DMA_M2M = DMA_M2M_Disable;



   /* Initialise the DMA */

   DMA_Init(DMA1_Channel1, &DMA_InitStructure);



   /* Enable DMA1_Channel3 Transfer Complete interrupt */

   DMA_ITConfig(DMA1_Channel1, DMA_IT_TC, ENABLE);



   /* Enable DMA1 Channel3_Tx */

   DMA_Cmd(DMA1_Channel1, ENABLE);

 }



 /**

   * @brief  Configures the Timers.

   * @param  wavePeriod (period of timer), preScaler (prescaler for timer)

   * @retval : None

   */

 void Timer_Configuration(uint16_t wavPeriod, uint16_t preScaler)

 {

 TIM_TimeBaseInitTypeDef  TIM_TimeBaseStructure;

 TIM_TimeBaseInitTypeDef TIM_TimeBaseInitStruct;



 //Time Base Configurations

 TIM_TimeBaseInitStruct.TIM_Period = 2000;

 TIM_TimeBaseInitStruct.TIM_Prescaler = 24000 - 1;

 TIM_TimeBaseInitStruct.TIM_ClockDivision = 0;

 TIM_TimeBaseInitStruct.TIM_CounterMode = TIM_CounterMode_Up;



 /* Time base configuration */

 TIM_TimeBaseStructure.TIM_Period = wavPeriod-1;

 TIM_TimeBaseStructure.TIM_Prescaler = preScaler-1;

 TIM_TimeBaseStructure.TIM_ClockDivision = 0;

 TIM_TimeBaseStructure.TIM_CounterMode = TIM_CounterMode_Up;



 /* Initialise the Timer */

 TIM_TimeBaseInit(TIM6, &TIM_TimeBaseStructure);

 TIM_TimeBaseInit(TIM2, &TIM_TimeBaseInitStruct);



 /* Need this to Trigger the DAC*/

 TIM_SelectOutputTrigger(TIM6, TIM_TRGOSource_Update);



 /*Enable TIMER2 interrupt*/

 TIM_ITConfig(TIM2, TIM_IT_Update, ENABLE);



 /* TIM6 enable counter */

 TIM_Cmd(TIM6, ENABLE);

 TIM_Cmd(TIM2, ENABLE);
 }



 /**

   * @brief  Configures the DAC

   * @param  None

   * @retval : None

   */


 void DAC_Configuration(void)

 {

 DAC_InitTypeDef DAC_InitStructure;



 /* Setup DAC */

 /* Default DAC struct init */

 DAC_StructInit(&DAC_InitStructure); /* Initiliase the struct */



 /* Change struct params */

 DAC_InitStructure.DAC_Trigger = DAC_Trigger_T6_TRGO;  /* Trigger on timer 6 */
 TIM6->ARR = G4Period-1;

 DAC_InitStructure.DAC_OutputBuffer = DAC_OutputBuffer_Enable; /* Enable
 output buffer -- can drive speaker */



 /* Init DAC */

 DAC_Init(DAC_Channel_1,&DAC_InitStructure); /* Initialise DAC */


 DAC_DMACmd(DAC_Channel_1, ENABLE); /* Enable DMA Request */

 DAC_Cmd(DAC_Channel_1, ENABLE); /* Enable the DAC */



 }
